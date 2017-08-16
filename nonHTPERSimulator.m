function packetErrorRate = nonHTPERSimulator(cfgNHT, chan, cfgChan, snr, ...
    maxNumErrors, maxNumPackets, enableFE, cfgLINC)
%nonHTPERSimulator Featured example helper function
%
%   Simulates the Non-HT transmit-receive link over a fading channel.

%   Copyright 2015-2016 The MathWorks, Inc.

% Waveform generation parameters
idleTime = 0;
numPkts = 1;
winTransTime = 1.6e-6; % No windowing

fs = helperSampleRate(cfgNHT);     % Baseband sampling rate

% Indices for accessing each field within the time-domain packet
ind = wlanFieldIndices(cfgNHT);

% Get the number of occupied subcarriers and FFT length
ofdmInfo = wlan.internal.wlanGetOFDMConfig(cfgNHT.ChannelBandwidth, 'Long', 'Legacy');
Nst = numel(ofdmInfo.DataIndices)+numel(ofdmInfo.PilotIndices); % Number of occupied subcarriers

% Create an instance of the AWGN channel per SNR point simulated
awgnChannel = comm.AWGNChannel;
awgnChannel.NoiseMethod = 'Signal to noise ratio (SNR)';
awgnChannel.SignalPower = 1;              % Unit power
awgnChannel.SNR = snr-10*log10(ofdmInfo.FFTLength/Nst); % Account for energy in nulls
awgnChannel.RandomStream = 'mt19937ar with seed'; % Local stream for repeatability 
release(awgnChannel);
reset(awgnChannel);


chDelay = 100; % arbitrary delay to account for all channel profiles

% Loop to simulate multiple packets

if cfgLINC.aFiltFlag == 1 && cfgLINC.flg  
    numPktErrStnd = 0;
    numPktErrLINC = 0;
    numPktErrLwTF = 0;
    numPkt = 0;
elseif cfgLINC.flg
    numPktErrStnd = 0;
    numPktErrLINC = 0;
    numPkt = 0;
else
    numPacketErrors = 0;
    numPkt = 0; % Index of packet transmitted
end

if cfgLINC.aFiltFlag == 1 && cfgLINC.flg
    while (numPktErrStnd<=maxNumErrors || numPktErrLINC<=maxNumErrors || numPktErrLwTF<=maxNumErrors) && numPkt<=maxNumPackets
        % Generate a packet waveform
        inpPSDU = randi([0 1], cfgNHT.PSDULength*8, 1); % PSDULength in bytes

        tx = wlanWaveformGenerator(inpPSDU,cfgNHT, 'IdleTime', idleTime,...
            'NumPackets', numPkts, 'WindowTransitionTime', winTransTime);

        % In the next  segment the idea is to have the same data in all paths no
        % matter if it is with or without LINC and it's transfer funciton
        % Make 2 paths for LINC and w/o LINC
        
        % Add trailing zeros to allow for channel delay
        padTx = [tx; zeros(chDelay, 1)];
        
        cfgLINC.aFiltFlag = 0;
        t = linc(padTx,cfgLINC);
        cfgLINC.aFiltFlag = 1;
        y = linc(padTx,cfgLINC); 
        maxDelay = max(t.delay,y.delay);

        % Add trailing zeros to allow for linc delay
        padTxStnd = [zeros(maxDelay,1);padTx];
        lengthStnd = length(padTxStnd);
        
        padTxLINC = [zeros(lengthStnd-length(t.do),1);t.do];
        padTxLwTF = [zeros(lengthStnd-length(y.do),1);y.do];
        if length(padTxLINC) == length(padTxLwTF)
            disp('match');
        end



        % Pass through HiperLAN/2 fading channel model
        if cfgChan.mPath
            rxStnd = filter(chan, padTxStnd);       % Reset channel to create different realizations
            %     disp(chan)
            reset(chan);
            %     disp(chan)
            rxLINC = filter(chan,padTxLINC);
            %     disp(chan)
            reset(chan);                        % Reset channel to create different realizations
            %     disp(chan)
            rxLwTF = filter(chan,padTxLwTF);
            %     disp(chan)
            reset(chan);                        % Reset channel to create different realizations
            %     disp(chan)
        else
            rxStnd = padTxStnd;
            rxLINC = padTxLINC;
            rxLwTF = padTxLwTF;
        end
        

        % Add white noise
        if cfgChan.awgn
            rxStnd = awgnChannel(rxStnd);
            %disp(awgnChannel);
            reset(awgnChannel);
            %disp(awgnChannel);

            rxLINC = awgnChannel(rxLINC);
            %disp(awgnChannel);
            reset(awgnChannel);
            %disp(awgnChannel);

            rxLwTF = awgnChannel(rxLwTF);
            %disp(awgnChannel);
            reset(awgnChannel);
            %disp(awgnChannel);
        end
        

        
        doStnd = true;
        doLINC = true;
        doLwTF = true;
        if enableFE
            % Packet detect
            pktOffsetStnd = wlanPacketDetect(rxStnd, cfgNHT.ChannelBandwidth);
            pktOffsetLINC = wlanPacketDetect(rxLINC, cfgNHT.ChannelBandwidth);
            pktOffsetLwTF = wlanPacketDetect(rxLwTF, cfgNHT.ChannelBandwidth);
            
            if isempty(pktOffsetStnd) % If empty no L-STF detected; packet error
                numPktErrStnd = numPktErrStnd+1;
                doStnd = false;                
            end
            
            if isempty(pktOffsetLINC)
               numPktErrLINC = numPktErrLINC+1;
               doLINC = false;
            end
            
            if isempty(pktOffsetLwTF)
               numPktErrLwTF = numPktErrLwTF+1;
               doLwTF = false;
            end
        
            if ~(doStnd || doLINC || doLwTF)
                numPkt = numPkt+1;
                continue; % Go to next loop iteration
            end
            
            if doStnd
                % Extract L-STF and perform coarse frequency offset correction
                lstfStnd = rxStnd(pktOffsetStnd+(ind.LSTF(1):ind.LSTF(2)),:);
                coarseFreqOffStnd = wlanCoarseCFOEstimate(lstfStnd, cfgNHT.ChannelBandwidth);
                rxStnd = helperFrequencyOffset(rxStnd, fs, -coarseFreqOffStnd);

                % Extract the Non-HT fields and determine start of L-LTF
                nonhtfieldsStnd = rxStnd(pktOffsetStnd+(ind.LSTF(1):ind.LSIG(2)),:);
                lltfIdxStnd = helperSymbolTiming(nonhtfieldsStnd, cfgNHT.ChannelBandwidth);
                
                % Synchronize the received waveform given the offset between the
                % expected start of the L-LTF and actual start of L-LTF
                pktOffsetStnd = pktOffsetStnd+lltfIdxStnd-double(ind.LLTF(1));
                
                % If no L-LTF detected or if packet detected outside the range of
                % expected delays from the channel modeling; packet error
                if isempty(lltfIdxStnd) || pktOffsetStnd<0 || pktOffsetStnd>chDelay
                    numPktErrStnd = numPktErrStnd+1;
                    doStnd = false;
                end
            end
            
            if doLINC
                % Extract L-STF and perform coarse frequency offset correction
                lstfLINC = rxLINC(pktOffsetLINC+(ind.LSTF(1):ind.LSTF(2)),:);
                coarseFreqOffLINC = wlanCoarseCFOEstimate(lstfLINC, cfgNHT.ChannelBandwidth);
                rxLINC = helperFrequencyOffset(rxLINC, fs, -coarseFreqOffLINC);

                % Extract the Non-HT fields and determine start of L-LTF
                nonhtfieldsLINC = rxLINC(pktOffsetLINC+(ind.LSTF(1):ind.LSIG(2)),:);
                lltfIdxLINC = helperSymbolTiming(nonhtfieldsLINC, cfgNHT.ChannelBandwidth);

                % Synchronize the received waveform given the offset between the
                % expected start of the L-LTF and actual start of L-LTF
                pktOffsetLINC = pktOffsetLINC+lltfIdxLINC-double(ind.LLTF(1));

                % If no L-LTF detected or if packet detected outside the range of
                % expected delays from the channel modeling; packet error            
                if isempty(lltfIdxLINC) || pktOffsetLINC<0 || pktOffsetLINC>chDelay
                    numPktErrLINC = numPktErrLINC+1;
                    doLINC = false;
                end
            end
            
            if doLwTF
                % Extract L-STF and perform coarse frequency offset correction
                lstfLwTF = rxLwTF(pktOffsetLwTF+(ind.LSTF(1):ind.LSTF(2)),:);
                coarseFreqOffLwTF = wlanCoarseCFOEstimate(lstfLwTF, cfgNHT.ChannelBandwidth);
                rxLwTF = helperFrequencyOffset(rxLwTF, fs, -coarseFreqOffLwTF);

                % Extract the Non-HT fields and determine start of L-LTF
                nonhtfieldsLwTF = rxLwTF(pktOffsetLwTF+(ind.LSTF(1):ind.LSIG(2)),:);
                lltfIdxLwTF = helperSymbolTiming(nonhtfieldsLwTF, cfgNHT.ChannelBandwidth);

                % Synchronize the received waveform given the offset between the
                % expected start of the L-LTF and actual start of L-LTF
                pktOffsetLwTF = pktOffsetLwTF+lltfIdxLwTF-double(ind.LLTF(1));

                % If no L-LTF detected or if packet detected outside the range of
                % expected delays from the channel modeling; packet error            
                if isempty(lltfIdxLwTF) || pktOffsetLwTF<0 || pktOffsetLwTF>chDelay
                    numPktErrLwTF = numPktErrLwTF+1;
                    doLwTF = false;
                end
            end
            
            if ~(doStnd || doLINC || doLwTF)
                numPkt = numPkt+1;
                continue; % Go to next loop iteration
            end
            
            if doStnd
                rxStnd = rxStnd(1+pktOffsetStnd:end,:);
                % Extract L-LTF and perform fine frequency offset correction
                lltfStnd = rxStnd(ind.LLTF(1):ind.LLTF(2),:);
                fineFreqOffStnd = wlanFineCFOEstimate(lltfStnd, cfgNHT.ChannelBandwidth);
                rxStnd = helperFrequencyOffset(rxStnd, fs, -fineFreqOffStnd);
            end
            if doLINC
            	rxLINC = rxLINC(1+pktOffsetLINC:end,:);
                % Extract L-LTF and perform fine frequency offset correction
                lltfLINC = rxLINC(ind.LLTF(1):ind.LLTF(2),:);
                fineFreqOffLINC = wlanFineCFOEstimate(lltfLINC, cfgNHT.ChannelBandwidth);
                rxLINC = helperFrequencyOffset(rxLINC, fs, -fineFreqOffLINC);
            end
            if doLwTF
            	rxLwTF = rxLwTF(1+pktOffsetLwTF:end,:);
                % Extract L-LTF and perform fine frequency offset correction
                lltfLwTF = rxLwTF(ind.LLTF(1):ind.LLTF(2),:);
                fineFreqOffLwTF = wlanFineCFOEstimate(lltfLwTF, cfgNHT.ChannelBandwidth);
                rxLwTF = helperFrequencyOffset(rxLwTF, fs, -fineFreqOffLwTF);
            end
        else
            
            % Directly offset the simulated channel filter delay
            chDelay = chan.ChannelFilterDelay; 
            chDelay = chDelay + t.delay;
            rxStnd = rxStnd(chDelay+1:end,:);
            rxLINC = rxLINC(chDelay+1:end,:);
            rxLwTF = rxLwTF(chDelay+1:end,:);
        end
        
        if doStnd
            % Extract L-LTF samples from the waveform, demodulate and perform
            % channel estimation
            lltfStnd = rxStnd(ind.LLTF(1):ind.LLTF(2),:);
            lltfDemodStnd = wlanLLTFDemodulate(lltfStnd, cfgNHT, 1);
            chanEstStnd = wlanLLTFChannelEstimate(lltfDemodStnd, cfgNHT);

            % Get estimate of the noise power from L-LTF
            nVarStnd = helperNoiseEstimate(lltfDemodStnd);

            % Extract Non-HT Data samples from the waveform and recover the PSDU
            nhtdataStnd = rxStnd(ind.NonHTData(1):ind.NonHTData(2),:);
            rxPSDUStnd = wlanNonHTDataRecover(nhtdataStnd, chanEstStnd, nVarStnd, cfgNHT);

            % Determine if any bits are in error, i.e. a packet error
            packetErrorStnd = any(biterr(inpPSDU, rxPSDUStnd));
            numPktErrStnd = numPktErrStnd+packetErrorStnd;
               
        end
        
        if doLINC
            % Extract L-LTF samples from the waveform, demodulate and perform
            % channel estimation
            lltfLINC = rxLINC(ind.LLTF(1):ind.LLTF(2),:);
            lltfDemodLINC = wlanLLTFDemodulate(lltfLINC, cfgNHT, 1);
            chanEstLINC = wlanLLTFChannelEstimate(lltfDemodLINC, cfgNHT);

            % Get estimate of the noise power from L-LTF
            nVarLINC = helperNoiseEstimate(lltfDemodLINC);

            % Extract Non-HT Data samples from the waveform and recover the PSDU
            nhtdataLINC = rxLINC(ind.NonHTData(1):ind.NonHTData(2),:);
            rxPSDULINC = wlanNonHTDataRecover(nhtdataLINC, chanEstLINC, nVarLINC, cfgNHT);

            % Determine if any bits are in error, i.e. a packet error
            packetErrorLINC = any(biterr(inpPSDU, rxPSDULINC));
            numPktErrLINC = numPktErrLINC+packetErrorLINC;   
        end
        if doLwTF
            % Extract L-LTF samples from the waveform, demodulate and perform
            % channel estimation
            lltfLwTF = rxLwTF(ind.LLTF(1):ind.LLTF(2),:);
            lltfDemodLwTF = wlanLLTFDemodulate(lltfLwTF, cfgNHT, 1);
            chanEstLwTF = wlanLLTFChannelEstimate(lltfDemodLwTF, cfgNHT);

            % Get estimate of the noise power from L-LTF
            nVarLwTF = helperNoiseEstimate(lltfDemodLwTF);

            % Extract Non-HT Data samples from the waveform and recover the PSDU
            nhtdataLwTF = rxLwTF(ind.NonHTData(1):ind.NonHTData(2),:);
            rxPSDULwTF = wlanNonHTDataRecover(nhtdataLwTF, chanEstLwTF, nVarLwTF, cfgNHT);

            % Determine if any bits are in error, i.e. a packet error
            packetErrorLwTF = any(biterr(inpPSDU, rxPSDULwTF));
            numPktErrLwTF = numPktErrLwTF+packetErrorLwTF;   
        end
        numPkt = numPkt+1; 
    end  
elseif cfgLINC.flg
    while (numPktErrStnd<=maxNumErrors || numPktErrLINC<=maxNumErrors) && numPkt<=maxNumPackets
        % Generate a packet waveform
        inpPSDU = randi([0 1], cfgNHT.PSDULength*8, 1); % PSDULength in bytes

        tx = wlanWaveformGenerator(inpPSDU,cfgNHT, 'IdleTime', idleTime,...
            'NumPackets', numPkts, 'WindowTransitionTime', winTransTime);

        % In the next  segment the idea is to have the same data in all paths no
        % matter if it is with or without LINC and it's transfer funciton
        % Make 2 paths for LINC and w/o LINC
        
        % Add trailing zeros to allow for channel delay
        padTx = [tx; zeros(chDelay, 1)];
        

        t = linc(padTx,cfgLINC);
        
        padTxLINC = t.do;
        padTxStnd = [zeros(t.delay,1); padTx; zeros(length(padTxLINC)-length(padTx)-t.delay,1)];
%         ratio = mean(abs(padTxStnd))/mean(abs(padTxLINC));
%         padTxLINC = padTxLINC*ratio;
%         
        % Add trailing zeros to allow for linc delay
        
        % Pass through HiperLAN/2 fading channel model
        if cfgChan.mPath
            rxStnd = filter(chan, padTxStnd);       % Reset channel to create different realizations
            %     disp(chan)
            reset(chan);
            %     disp(chan)
            rxLINC = filter(chan,padTxLINC);
            %     disp(chan)
            reset(chan);                        % Reset channel to create different realizations
            %     disp(chan)           
        else
            rxStnd = padTxStnd;
            rxLINC = padTxLINC;
        end

        % Add noise
        if cfgChan.awgn
            rxStnd = awgnChannel(rxStnd);
            %disp(awgnChannel);
            reset(awgnChannel);
            %disp(awgnChannel);

            rxLINC = awgnChannel(rxLINC);
            %disp(awgnChannel);
            reset(awgnChannel);
            %disp(awgnChannel);    
        end
        

        doStnd = true;
        doLINC = true;
        if enableFE
            % Packet detect
            pktOffsetStnd = wlanPacketDetect(rxStnd, cfgNHT.ChannelBandwidth);
            pktOffsetLINC = wlanPacketDetect(rxLINC, cfgNHT.ChannelBandwidth);
            
            if isempty(pktOffsetStnd) % If empty no L-STF detected; packet error
                numPktErrStnd = numPktErrStnd+1;
                doStnd = false;                
            end
            
            if isempty(pktOffsetLINC)
               numPktErrLINC = numPktErrLINC+1;
               doLINC = false;
            end
        
            if ~(doStnd || doLINC)
                numPkt = numPkt+1;
                continue; % Go to next loop iteration
            end
            
            if doStnd
                % Extract L-STF and perform coarse frequency offset correction
                lstfStnd = rxStnd(pktOffsetStnd+(ind.LSTF(1):ind.LSTF(2)),:);
                coarseFreqOffStnd = wlanCoarseCFOEstimate(lstfStnd, cfgNHT.ChannelBandwidth);
                rxStnd = helperFrequencyOffset(rxStnd, fs, -coarseFreqOffStnd);

                % Extract the Non-HT fields and determine start of L-LTF
                nonhtfieldsStnd = rxStnd(pktOffsetStnd+(ind.LSTF(1):ind.LSIG(2)),:);
                lltfIdxStnd = helperSymbolTiming(nonhtfieldsStnd, cfgNHT.ChannelBandwidth);
                
                % Synchronize the received waveform given the offset between the
                % expected start of the L-LTF and actual start of L-LTF
                pktOffsetStnd = pktOffsetStnd+lltfIdxStnd-double(ind.LLTF(1));
                
                % If no L-LTF detected or if packet detected outside the range of
                % expected delays from the channel modeling; packet error
                if isempty(lltfIdxStnd) || pktOffsetStnd<0 || pktOffsetStnd>chDelay
                    numPktErrStnd = numPktErrStnd+1;
                    doStnd = false;
                end
            end
            
            if doLINC
                % Extract L-STF and perform coarse frequency offset correction
                lstfLINC = rxLINC(pktOffsetLINC+(ind.LSTF(1):ind.LSTF(2)),:);
                coarseFreqOffLINC = wlanCoarseCFOEstimate(lstfLINC, cfgNHT.ChannelBandwidth);
                rxLINC = helperFrequencyOffset(rxLINC, fs, -coarseFreqOffLINC);

                % Extract the Non-HT fields and determine start of L-LTF
                nonhtfieldsLINC = rxLINC(pktOffsetLINC+(ind.LSTF(1):ind.LSIG(2)),:);
                lltfIdxLINC = helperSymbolTiming(nonhtfieldsLINC, cfgNHT.ChannelBandwidth);

                % Synchronize the received waveform given the offset between the
                % expected start of the L-LTF and actual start of L-LTF
                pktOffsetLINC = pktOffsetLINC+lltfIdxLINC-double(ind.LLTF(1));

                % If no L-LTF detected or if packet detected outside the range of
                % expected delays from the channel modeling; packet error            
                if isempty(lltfIdxLINC) || pktOffsetLINC<0 || pktOffsetLINC>chDelay
                    numPktErrLINC = numPktErrLINC+1;
                    doLINC = false;
                end
            end
            
            if ~(doStnd || doLINC)
                numPkt = numPkt+1;
                continue; % Go to next loop iteration
            end
            
            if doStnd
                rxStnd = rxStnd(1+pktOffsetStnd:end,:);
                % Extract L-LTF and perform fine frequency offset correction
                lltfStnd = rxStnd(ind.LLTF(1):ind.LLTF(2),:);
                fineFreqOffStnd = wlanFineCFOEstimate(lltfStnd, cfgNHT.ChannelBandwidth);
                rxStnd = helperFrequencyOffset(rxStnd, fs, -fineFreqOffStnd);
            end
            if doLINC
            	rxLINC = rxLINC(1+pktOffsetLINC:end,:);
                % Extract L-LTF and perform fine frequency offset correction
                lltfLINC = rxLINC(ind.LLTF(1):ind.LLTF(2),:);
                fineFreqOffLINC = wlanFineCFOEstimate(lltfLINC, cfgNHT.ChannelBandwidth);
                rxLINC = helperFrequencyOffset(rxLINC, fs, -fineFreqOffLINC);
            end
        else
            
            % Directly offset the simulated channel filter delay
            chDelay = chan.ChannelFilterDelay; 
            chDelay = chDelay + t.delay;
            rxStnd = rxStnd(chDelay+1:end,:);
            rxLINC = rxLINC(chDelay+1:end,:);
        end
        
        if doStnd
            % Extract L-LTF samples from the waveform, demodulate and perform
            % channel estimation
            lltfStnd = rxStnd(ind.LLTF(1):ind.LLTF(2),:);
            lltfDemodStnd = wlanLLTFDemodulate(lltfStnd, cfgNHT, 1);
            chanEstStnd = wlanLLTFChannelEstimate(lltfDemodStnd, cfgNHT);

            % Get estimate of the noise power from L-LTF
            nVarStnd = helperNoiseEstimate(lltfDemodStnd);

            % Extract Non-HT Data samples from the waveform and recover the PSDU
            nhtdataStnd = rxStnd(ind.NonHTData(1):ind.NonHTData(2),:);
            rxPSDUStnd = wlanNonHTDataRecover(nhtdataStnd, chanEstStnd, nVarStnd, cfgNHT);

            % Determine if any bits are in error, i.e. a packet error
            packetErrorStnd = any(biterr(inpPSDU, rxPSDUStnd));
            numPktErrStnd = numPktErrStnd+packetErrorStnd;
               
        end
        
        if doLINC
            % Extract L-LTF samples from the waveform, demodulate and perform
            % channel estimation
            lltfLINC = rxLINC(ind.LLTF(1):ind.LLTF(2),:);
            lltfDemodLINC = wlanLLTFDemodulate(lltfLINC, cfgNHT, 1);
            chanEstLINC = wlanLLTFChannelEstimate(lltfDemodLINC, cfgNHT);

            % Get estimate of the noise power from L-LTF
            nVarLINC = helperNoiseEstimate(lltfDemodLINC);

            % Extract Non-HT Data samples from the waveform and recover the PSDU
            nhtdataLINC = rxLINC(ind.NonHTData(1):ind.NonHTData(2),:);
            rxPSDULINC = wlanNonHTDataRecover(nhtdataLINC, chanEstLINC, nVarLINC, cfgNHT);

            % Determine if any bits are in error, i.e. a packet error
            packetErrorLINC = any(biterr(inpPSDU, rxPSDULINC));
            numPktErrLINC = numPktErrLINC+packetErrorLINC;   
        end
        numPkt = numPkt+1; 
    end  
else
    while numPacketErrors<=maxNumErrors && numPkt<=maxNumPackets
        % Generate a packet waveform
        inpPSDU = randi([0 1], cfgNHT.PSDULength*8, 1); % PSDULength in bytes

        tx = wlanWaveformGenerator(inpPSDU,cfgNHT, 'IdleTime', idleTime,...
            'NumPackets', numPkts, 'WindowTransitionTime', winTransTime);

        % Add trailing zeros to allow for channel delay
        padTx = [tx; zeros(chDelay, 1)];
        
        % Pass through HiperLAN/2 fading channel model
        if cfgChan.mPath
            rx = filter(chan, padTx);
            reset(chan);    % Reset channel to create different realizations
        else
            rx = padTx;
        end
        % Add noise
        if cfgChan.awgn
            rx = awgnChannel(rx);
        end
        
        if enableFE
            % Packet detect
            pktOffset = wlanPacketDetect(rx, cfgNHT.ChannelBandwidth);
            if isempty(pktOffset) % If empty no L-STF detected; packet error
                numPacketErrors = numPacketErrors+1;
                numPkt = numPkt+1;
                continue; % Go to next loop iteration
            end

            % Extract L-STF and perform coarse frequency offset correction
            lstf = rx(pktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
            coarseFreqOff = wlanCoarseCFOEstimate(lstf, cfgNHT.ChannelBandwidth);
            rx = helperFrequencyOffset(rx, fs, -coarseFreqOff);

            % Extract the Non-HT fields and determine start of L-LTF
            nonhtfields = rx(pktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
            lltfIdx = helperSymbolTiming(nonhtfields, cfgNHT.ChannelBandwidth);

            % Synchronize the received waveform given the offset between the
            % expected start of the L-LTF and actual start of L-LTF
            pktOffset = pktOffset+lltfIdx-double(ind.LLTF(1));
            % If no L-LTF detected or if packet detected outside the range of
            % expected delays from the channel modeling; packet error
            if isempty(lltfIdx) || pktOffset<0 || pktOffset>chDelay
                numPacketErrors = numPacketErrors+1;
                numPkt = numPkt+1;
                continue; % Go to next loop iteration
            end
            rx = rx(1+pktOffset:end,:);

            % Extract L-LTF and perform fine frequency offset correction
            lltf = rx(ind.LLTF(1):ind.LLTF(2),:);
            fineFreqOff = wlanFineCFOEstimate(lltf, cfgNHT.ChannelBandwidth);
            rx = helperFrequencyOffset(rx, fs, -fineFreqOff);
        else
            % Directly offset the simulated channel filter delay
            chDelay = chan.ChannelFilterDelay;
            rx = rx(chDelay+1:end,:);

        end

        % Extract L-LTF samples from the waveform, demodulate and perform
        % channel estimation
        lltf = rx(ind.LLTF(1):ind.LLTF(2),:);
        lltfDemod = wlanLLTFDemodulate(lltf, cfgNHT, 1);
        chanEst = wlanLLTFChannelEstimate(lltfDemod, cfgNHT);

        % Get estimate of the noise power from L-LTF
        nVar = helperNoiseEstimate(lltfDemod);

        % Extract Non-HT Data samples from the waveform and recover the PSDU
        nhtdata = rx(ind.NonHTData(1):ind.NonHTData(2),:);
        rxPSDU = wlanNonHTDataRecover(nhtdata, chanEst, nVar, cfgNHT);

        % Determine if any bits are in error, i.e. a packet error
        packetError = any(biterr(inpPSDU, rxPSDU));
        numPacketErrors = numPacketErrors+packetError;
        numPkt = numPkt+1;
    end
end

% Calculate packet error rate (PER) at SNR point
if cfgLINC.aFiltFlag == 1 && cfgLINC.flg
    packetErrorRate.Stnd = numPktErrStnd/numPkt;
    disp(['CBW' cfgNHT.ChannelBandwidth(4:end) ', SNR ' num2str(snr) ...
    ' completed after ' num2str(numPkt) ' packets, PER: ' ...
    num2str(packetErrorRate.Stnd) ', LINC is OFF']);

    packetErrorRate.LINC = numPktErrLINC/numPkt;
    disp(['CBW' cfgNHT.ChannelBandwidth(4:end) ', SNR ' num2str(snr) ...
    ' completed after ' num2str(numPkt) ' packets, PER: ' ...
    num2str(packetErrorRate.LINC) ', LINC is  ON']);

    packetErrorRate.LwTF = numPktErrLwTF/numPkt;
    disp(['CBW' cfgNHT.ChannelBandwidth(4:end) ', SNR ' num2str(snr) ...
    ' completed after ' num2str(numPkt) ' packets, PER: ' ...
    num2str(packetErrorRate.LwTF) ', LwTF is  ON']);
elseif cfgLINC.flg
    packetErrorRate.Stnd = numPktErrStnd/numPkt;
    disp(['CBW' cfgNHT.ChannelBandwidth(4:end) ', SNR ' num2str(snr) ...
    ' completed after ' num2str(numPkt) ' packets, PER: ' ...
    num2str(packetErrorRate.Stnd) ', LINC is OFF']);
    packetErrorRate.LINC = numPktErrLINC/numPkt;
    disp(['CBW' cfgNHT.ChannelBandwidth(4:end) ', SNR ' num2str(snr) ...
    ' completed after ' num2str(numPkt) ' packets, PER: ' ...
    num2str(packetErrorRate.LINC) ', LINC is  ON']);
else
    packetErrorRate.Stnd = numPacketErrors/numPkt;
    disp(['CBW' cfgNHT.ChannelBandwidth(4:end) ', SNR ' num2str(snr) ...
    ' completed after ' num2str(numPkt) ' packets, PER: ' ...
    num2str(packetErrorRate.Stnd) ', LINC is OFF']);
end


end