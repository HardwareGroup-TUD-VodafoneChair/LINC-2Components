function t = linc(OFDMburst,cfg)
% I/Q baseband simulator of combiner network
%
%  IQ -> Upsampling -> LP Filter -> Theta/Phi -> 
%     -> Phase Mod. -> Combiner Baseband Filters -> Combiner -> 
%     -> Downsamling -> Error Analysis
%
% Author: E.Matus, Stefan Damjancevic
% Last edit: 7/18/17 SD
% 

%--------------------%
% Configuration      %
%--------------------%
t.cfg.mode       = cfg.mode;               % outphasing mode for a(t)=A(t)cos(wc.t+phi(t)) : 
                                    %                   1:theta=arcsin(abs(A(t)); 
                                    %                   2;theta=arccos(abs(A(t)) 
                                    %                   3;theta1=arcsin(real(a(t)),theta2=arcsin(imag(a(t))
%t.cfg.src        = 7;               % 1:BPSK, 2:harmonic, 3:LP-noise, 4:Rectangle, 5:Ramp, 6:Duo-harmonic, 7:OFDM, 8.Kronecker
t.cfg.fs         = cfg.fs;           % IQ sampling frequency 20e6 or 10e6
t.cfg.L          = length(OFDMburst);           % IQ vector length
%t.cfg.cp         = cfg.cp;          % cyclic prefix
t.cfg.sps        = cfg.sps;          % IQ - oversampling; default was 4;
t.cfg.aFiltFlag  = cfg.aFiltFlag;    % use combiner transfer function
t.cfg.lpFiltFlag = cfg.lpFiltFlag;   % Spectrum shaping filter
%t.cfg.errType    = 2;               % 1:comb-rrcFilt, 2:irrcFilt-iq; 




t.cfgPlot.iq        = 0;
t.cfgPlot.rrcFilt   = 0;
t.cfgPlot.lpFilt    = 0;
t.cfgPlot.angle     = 0;
t.cfgPlot.pm        = 0;
t.cfgPlot.aFilt     = 0;
t.cfgPlot.dFilt     = 0;
t.cfgPlot.comb      = 0;
t.cfgPlot.irrcFilt  = 0;
t.cfgPlot.err       = 0;

%--------------------%
% Initialization     % 
%--------------------%
figId       = 1;
t.cfg.N     = t.cfg.L * t.cfg.sps;        % length of oversampled symbol


%--------------------%
% Source IQ sapmples %
%--------------------%
t.iq.name        = 'IQ samples';
t.iq.plot        = t.cfgPlot.iq;
t.iq.figId       = figId; figId = figId+1;

t.iq.L           = t.cfg.L;                          % vector length of IQ signal

t.iq.do = OFDMburst;

% t.iq.max = max(abs(t.iq.do));
% t.iq.do = 0.9 * t.iq.do/t.iq.max;          % scale amplitude |iq|<0.9


%----------------------------------------------------------------------%
% Upsampling using root raised cosine filter                           %
% Ref: P.Xio - Transmit and Receive Filter Design For OFDM Based WLAN  %
% Systems.                                                             %
%----------------------------------------------------------------------%
t.rrcFilt.name    = 'Upsampling RRC Filter';
t.rrcFilt.plot    = t.cfgPlot.rrcFilt;
t.rrcFilt.figId   = figId; figId = figId+1;

t.rrcFilt.rolloff = 0.2;
t.rrcFilt.span    = 6;
t.rrcFilt.sps     = t.cfg.sps;
t.rrcFilt.b       = rcosdesign(t.rrcFilt.rolloff, t.rrcFilt.span, t.rrcFilt.sps,'normal');
% t.rrcFilt.b       = t.rrcFilt.b / sum(t.rrcFilt.b);                                             % filter scaling
t.rrcFilt.delay   = round(mean(grpdelay(t.rrcFilt.b)));                                         % mean group delay
t.rrcFilt.L       = t.iq.L * t.rrcFilt.sps;
% t.rrcFilt.n1      = (t.iq.n1-1) * t.rrcFilt.sps + 1 + length(t.rrcFilt.b) - 1; % t.rrcFilt.delay;                  
% t.rrcFilt.n2      = t.rrcFilt.n1 + t.rrcFilt.L - 1;

t.rrcFilt.do      = upfirdn(real(t.iq.do), t.rrcFilt.b, t.rrcFilt.sps, 1);                          % oversampled real
t.rrcFilt.do      = complex(t.rrcFilt.do, upfirdn(imag(t.iq.do), t.rrcFilt.b, t.rrcFilt.sps, 1));   % oversampled imag
t.rrcFilt.do      = t.rrcFilt.do * t.rrcFilt.sps;                                                   % scales amplitude


% tmp             = max(min(real(t.rrcFilt.do),1),-1);                                              % clipping to <-1,1>
% t.rrcFilt.do    = complex(tmp, max(min(imag(t.rrcFilt.do),1),-1));                                % clipping to <-1,1>

%------------------------------------------------%
% LP FIR Filter                                  %
% Simulate spectral mask                         %
% 802.11 Spectral mask (symmetry to 0 MHz)       %
% 20MHz channel:                                 %
% MHz:      0, 9,  11,  20,  30                  %
% PSD [dB]: 0, 0, -20, -28, -45                  %
%                                                %
% 40MHz channel:                                 %
% MHz:      0, 19,  21,  40,  60                 %
% PSD [dB]: 0,  0, -20, -28, -45                 %
%                                                %
% Spectral flatness:                             %
% SC 1 to 16 (-1 to -16): <2dB to avg SC 1-16    %
% SC 17 to 28 (-17 to -28): +2/-4dB to SC 1-16   %
%------------------------------------------------%
t.lpFilt.name   = 'IQ Spectral mask';
t.lpFilt.plot   = t.cfgPlot.lpFilt;
t.lpFilt.figId  = figId; figId = figId+1;

t.lpFilt.sps    = t.rrcFilt.sps;
t.lpFilt.normFc = .85/t.lpFilt.sps;
t.lpFilt.order  = 32;
t.lpFilt.beta   = 6;
lpFilt  = firls(t.lpFilt.order, [0 t.lpFilt.normFc t.lpFilt.normFc 1],[1 1 0 0]);
lpFilt  = lpFilt .* kaiser(t.lpFilt.order+1, t.lpFilt.beta)';

if t.cfg.lpFiltFlag == 1
    t.lpFilt.b  = lpFilt;   % user defined FIR
else
    t.lpFilt.b  = [1];      % Bypass filter
end

t.lpFilt.b      = t.lpFilt.b/sum(t.lpFilt.b);           % normalize DC=1
t.lpFilt.delay  = round(mean(grpdelay(t.lpFilt.b)));    % group delay
t.lpFilt.L      = t.rrcFilt.L;
% t.lpFilt.n1     = t.rrcFilt.n1 + length(t.lpFilt.b) - 1;   
% t.lpFilt.n2     = t.lpFilt.n1  + t.lpFilt.L - 1;
t.lpFilt.do     = filter(t.lpFilt.b, 1, t.rrcFilt.do);

t.lpFilt.max    = max(abs(t.lpFilt.do));
t.lpFilt.do     = 0.96*t.lpFilt.do/t.lpFilt.max;        % scale amplitude |do|< 1

%---------------------------------%
% IQ to Theta/Phi transformation  %
%---------------------------------%
t.angle.name  = 'Theta/Phi';
t.angle.plot  = t.cfgPlot.angle;
t.angle.figId = figId; figId = figId+1;

t.angle.L       = t.lpFilt.L;
% t.angle.n1      = t.lpFilt.n1;
% t.angle.n2      = t.lpFilt.n2;
tmp             = t.lpFilt.do; % / max(abs(t.lpFilt.do));   % amplitude scaling to <-1,1>
t.angle.phi     = zeros(size(tmp));

if t.cfg.mode == 1
    t.angle.phi     = unwrap( angle(tmp));
    t.angle.theta	= asin(abs(tmp));
    t.angle.alpha   = t.angle.phi + t.angle.theta;
    t.angle.beta    = t.angle.phi - t.angle.theta;
elseif t.cfg.mode == 2
    t.angle.phi     = unwrap( angle(tmp));
    t.angle.theta   = acos(abs(tmp));
    t.angle.alpha   = t.angle.phi + t.angle.theta;
    t.angle.beta    = t.angle.phi - t.angle.theta;
elseif t.cfg.mode == 3
    t.angle.theta{1} = asin(real(tmp));
    t.angle.theta{2} = asin(imag(tmp));
end


%-----------------%
% Phase Modulator %
%-----------------%
t.pm.name = 'Phase Modulator';
t.pm.plot = t.cfgPlot.pm;
t.pm.figId = figId; figId =  figId+1;

t.pm.L  = t.angle.L;
% t.pm.n1 = t.angle.n1;
% t.pm.n2 = t.angle.n2;

%if(any(isnan(t.angle.alpha))|| any(isinf(t.angle.alpha)) || any(isnan(t.angle.beta))|| any(isinf(t.angle.beta)) || any(isnan(t.angle.theta))|| any(isinf(t.angle.theta)))
if any(~isfinite(t.angle.alpha)) || any(~isfinite(t.angle.beta)) || any(~isfinite(t.angle.theta)) % cannot detect the error that causes complex to stop
    disp('pause')
end
if t.cfg.mode == 1
    t.pm.do{1} =  0.5*complex( +sin(t.angle.alpha ),-cos(t.angle.alpha ));
    t.pm.do{2} =  0.5*complex( -sin(t.angle.beta ),  +cos(t.angle.beta ));
elseif t.cfg.mode == 2
    t.pm.do{1} =  0.5*complex( +cos(t.angle.alpha), +sin(t.angle.alpha));
    t.pm.do{2} =  0.5*complex( +cos(t.angle.beta), +sin(t.angle.beta)); % without complex it passes 
elseif t.cfg.mode == 3
    t.pm.do1{1} = 0.5*complex( +sin(t.angle.theta{1}), -cos(t.angle.theta{1}));
    t.pm.do1{2} = 0.5*complex( +sin(t.angle.theta{1}), +cos(t.angle.theta{1}));
    t.pm.do2{1} = 0.5*complex( -cos(t.angle.theta{2}), -sin(t.angle.theta{2}));
    t.pm.do2{2} = 0.5*complex( +cos(t.angle.theta{2}), -sin(t.angle.theta{2}));
    t.pm.do{1}  = t.pm.do1{1} - t.pm.do2{1};
    t.pm.do{2}  = t.pm.do1{2} - t.pm.do2{2};
end

%-------------------------%
% Analog combiner filters %
%-------------------------%
t.aFilt.name   = 'Analog Filter';
t.aFilt.plot   = t.cfgPlot.aFilt;
t.aFilt.figId  = figId; figId = figId+1;

t.aFilt.wc     = 2*pi*4.3e9;     % filter central frequency
t.aFilt.bw     = 2*pi*200e6;     % Band width
t.aFilt.ord    = 2;              % filter order
t.aFilt.Rp     = 5;              % ripple in passband [db]
t.aFilt.Rs     = 10;             % ripple in stopband [db]
t.aFilt.ai     = {1,1};          % initialize TF coeffs
t.aFilt.bi     = {1 1};          % initialize TF coeffs
t.aFilt.flag   = t.cfg.aFiltFlag;

if t.aFilt.flag == 1
    %[t.aFilt.bi{1},t.aFilt.ai{1}]  = ellip(t.aFilt.ord, t.aFilt.Rp, t.aFilt.Rs, t.aFilt.wc, 's');
    %[t.aFilt.bi{1},t.aFilt.ai{1}]  = butter(t.aFilt.ord, [t.aFilt.wc-t.aFilt.bw/2, t.aFilt.wc+t.aFilt.bw/2],'bandpass', 's');  
    %[t.aFilt.bi{1},t.aFilt.ai{1}]  = butter(t.aFilt.ord, t.aFilt.wc, 'low', 's');
    
    % test filters of 2-stage combiner
    %[t.aFilt.bi{1},t.aFilt.ai{1}]   = butter(t.aFilt.ord, t.aFilt.wc*1.01, 'low', 's');
    %[t.aFilt.bi{2},t.aFilt.ai{2}]   = butter(t.aFilt.ord, t.aFilt.wc, 'low', 's');
    
    % 2-stage combiner transfer function coefficients H1(s) and H2(s)
    t.aFilt.bi{1} = ...
       [ 0.000000000000000000000000000000051771998331461798199055311048463760804936998182155, ...
         0.0000000000000000000027981263722122299990077131579690065503229666230931, ...
         0.000000000037791277785615200347156690965873043416833887420125, ...
         0];
    t.aFilt.ai{1} = ...
       [ 0.00000000000000000000000000000010354399666292399051330371345172871890881929665322, ...
         0.000000000000000000008394379116636679464593756106210948251089526357958, ...
         0.00000000022520670317556500742491275528082191742518247679072, ...
         2.0];
    t.aFilt.bi{2} = ...
       [ 0, ...
         0.000000000000000000068497284777134801413518665844344755515983742071674, ...
         0.0000000037020796016551999620566426580205512952836954809754, ...
         50.0];
    t.aFilt.ai{2} = ...
       [ 0.0000000000000000000000000000049660961965417501755997792871589133147576336439725, ...
         0.0000000000000000004081906880377590195033748025231934467072055241311, ...
         0.000000011106238804965600299760234250575491321555432477908, ...
         100.0];

else
    t.aFilt.bi{1}   = [0 0 1 ];
    t.aFilt.ai{1}   = [0 0 1 ];
    t.aFilt.bi{2}   = [0 0 1 ];
    t.aFilt.ai{2}   = [0 0 1 ];
end


% Transfer functions
t.aFilt.tf0{1}   = tf(t.aFilt.bi{1}, t.aFilt.ai{1});
t.aFilt.tf0{2}   = tf(t.aFilt.bi{2}, t.aFilt.ai{2});
 
% Freq. shift of analog filter bi + b1.s + b2.s^2 +... by s0,
% b transformed to b' by subst s->(s+s0) according to:
%          k         k+1              k+2
% b(k)' = ( )b(k) + (   )b(k+1).s0 + (  )b(k+2).s0^2 ...
%          0          1                2
t.aFilt.w0   = -t.aFilt.wc;                   % freq. shift omega_0
t.aFilt.b{1} = zeros(size(t.aFilt.bi{1}));    % numerator freq. shifted filter
t.aFilt.a{1} = zeros(size(t.aFilt.ai{1}));    % denominator freq. shifted filter
t.aFilt.b{2} = zeros(size(t.aFilt.bi{2}));    % numerator freq. shifted filter
t.aFilt.a{2} = zeros(size(t.aFilt.ai{2}));    % denominator freq. shifted filter

oi = {'bi','ai'};                   % input coef field names
oo = {'b' ,'a' };                   % output coef field names
for m = 1:2,                        % for numerator b and denominator a
    in  = t.aFilt.(oi{m});  
    out = t.aFilt.(oo{m});  
    for n = 1:length(in),           % over all filters
        x = in{n};                  % coef vector
        N = length(x);          
        p = pascal(N+1);
        for k=1:N,
            v1 = x(1:end-k+1);
            v2 = fliplr(p(k,1:N-k+1));
            out{n}(end-k+1) = polyval(v1.*v2,(-j*t.aFilt.w0));
        end
    end
    t.aFilt.(oo{m})=out;
end

% Transfer function of freq. shifted filter
t.aFilt.tf{1} = tf(t.aFilt.b{1}, t.aFilt.a{1});
t.aFilt.tf{2} = tf(t.aFilt.b{2}, t.aFilt.a{2});

%-----------------------------------------%
% Digital representation of analog filter %
%-----------------------------------------%
t.bbFilt.name  = 'Digital Filter';
t.bbFilt.plot  = t.cfgPlot.dFilt;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
t.bbFilt.figId = figId; figId = figId+1;

t.bbFilt.fs    = t.aFilt.bw/2/pi * t.rrcFilt.sps;

if t.aFilt.flag ==1
    [t.bbFilt.b{1}, t.bbFilt.a{1}] = bilinear(t.aFilt.b{1}, t.aFilt.a{1}, t.bbFilt.fs); %, t.aFilt.wc/2/pi);
    [t.bbFilt.b{2}, t.bbFilt.a{2}] = bilinear(t.aFilt.b{2}, t.aFilt.a{2}, t.bbFilt.fs); %, t.aFilt.wc/2/pi);
else
    t.bbFilt.b{1} = [1,0,0];
    t.bbFilt.a{1} = [1,0,0];
    t.bbFilt.b{2} = [1,0,0];
    t.bbFilt.a{2} = [1,0,0];
end

g = grpdelay(t.bbFilt.b{1},t.bbFilt.a{1});
t.bbFilt.delay = ceil(mean(g(isfinite(g))));    % group delay, assume equal delay of all filters
t.bbFilt.L     = t.pm.L;
% t.bbFilt.n1    = t.pm.n1  + t.bbFilt.delay;
% t.bbFilt.n2    = t.bbFilt.n1 + t.bbFilt.L - 1;
if 1
    t.bbFilt.do{1} = filter(t.bbFilt.b{1}, t.bbFilt.a{1}, t.pm.do{1});     
    t.bbFilt.do{2} = filter(t.bbFilt.b{2}, t.bbFilt.a{2}, t.pm.do{2});     
else
    t.bbFilt.do{1} = filter(t.bbFilt.b{1}, t.bbFilt.a{1}, t.pm.do12{1});     
    t.bbFilt.do{2} = filter(t.bbFilt.b{2}, t.bbFilt.a{2}, t.pm.do12{2});     
end

%-----------%
% Combiner  %
%-----------%
t.comb.name = 'Combiner';
t.comb.plot = t.cfgPlot.comb;
t.comb.figId = figId; figId = figId+1;

t.comb.L  = t.bbFilt.L;
% t.comb.n1 = t.bbFilt.n1;
% t.comb.n2 = t.comb.n1 + t.comb.L - 1;
t.comb.do = t.bbFilt.do{1} + t.bbFilt.do{2};     

%--------------------------------%
% Downsampling using RRC filter  %
%--------------------------------%
t.irrcFilt.name = 'Downsampling inverse RRC filter';
t.irrcFilt.plot = t.cfgPlot.irrcFilt;
t.irrcFilt.figId= figId; figId = figId+1;

t.irrcFilt.sps  = t.rrcFilt.sps;                        % sample per symbol
t.irrcFilt.b    = t.rrcFilt.b;                          % filter coefs.
t.irrcFilt.delay= mean(grpdelay(t.irrcFilt.b));         % mean group delay
t.irrcFilt.L    = fix(t.comb.L/t.irrcFilt.sps);
% t.irrcFilt.n1   = t.comb.n1 + length(t.irrcFilt.b) - 1; % t.irrcFilt.delay;
% t.irrcFilt.n1   = ceil(t.irrcFilt.n1/t.irrcFilt.sps);
% t.irrcFilt.n2   = t.irrcFilt.n1 + t.irrcFilt.L - 1;

% downsample combiner output
t.irrcFilt.do = upfirdn(real(t.comb.do), t.irrcFilt.b, 1, t.irrcFilt.sps);                         % downsampled real
t.irrcFilt.do = complex(t.irrcFilt.do, upfirdn(imag(t.comb.do), t.irrcFilt.b, 1, t.irrcFilt.sps)); % downsampled imag


t.do = 1/(0.96/t.lpFilt.max) * t.irrcFilt.do;                  % scale back signal
t.delay = delay(t.do,OFDMburst);


%-----------------------%
% Error Analysis Block  %
%-----------------------%
%{
t.err.name  = 'Error Analysis';
t.err.plot  = t.cfgPlot.err;
t.err.figId = figId; figId = figId+1;
t.err.type  = t.cfg.errType;        % 1:comb-rrcFilt, 2:irrcFilt-iq; 
t.err.sps   = 1;                    % default

if t.err.type == 1
    t.err.di1   = t.rrcFilt.do;
    t.err.di1   = t.rrcFilt.cp;
    t.err.di2   = t.comb.do;
    t.err.n11   = t.rrcFilt.n1;
    t.err.n12   = t.comb.n1;
    t.err.L     = t.rrcFilt.L;
    t.err.sps   = t.rrcFilt.sps;
    t.err.idxSc = ones(1,t.err.L);
    t.err.ylabel= '|e|=|comb-rrcFilt|';
elseif t.err.type == 2
    t.err.di1   = t.iq.do;
    t.err.cp1   = t.iq.cp;
    t.err.di2   = t.irrcFilt.do;
    t.err.n11   = t.iq.n1;
    t.err.n12   = t.irrcFilt.n1;
    t.err.L     = t.iq.L;
    t.err.idxSc = t.iq.idxSc;
    t.err.ylabel= '|e(n)|=|y(n)-x(n)|';
end

t.err.idx1   = 1:length(t.err.di1); 
t.err.idx2   = 1:length(t.err.di2); 

% fit input and output sequences: shift and scale 
% TODO: implement better data fitting algorithm
shift = delay(t.err.di2(t.err.idx2), t.err.di1(t.err.idx1));
t.err.do2       = t.err.di2(t.err.n12:(t.err.n12+t.err.L-1));
t.err.do1       = t.err.di1( (t.err.n12:(t.err.n12+t.err.L-1)) - shift );
gain  = t.err.do1 ./ t.err.do2;
%idx   = find(isfinite(gain) & (abs(gain) > mean(abs(gain))/100  ) );
idx   = find(isfinite(gain) & (abs(gain) > 0) );
scale = mean(gain(idx));
 
% Equalizer
t.err.do2       = scale * t.err.do2;
%scale2          = sum(abs(t.err.do2).^2)/sum(abs(t.err.do1).^2);
%t.err.do2       = t.err.do2/sqrt(scale2);

t.err.eq1       = fft(t.err.do1)./fft(t.err.di1((1:t.err.L)+t.err.cp1));
t.err.idxSc0    = setdiff(1:t.err.L, t.err.idxSc);  % zero Subcarriers
t.err.eq1(t.err.idxSc0) = 1;                        % unused subcarriers

t.err.eq2       = fft(t.err.do2)./fft(t.err.di1((1:t.err.L)+t.err.cp1));
t.err.evm       = fft(t.err.do2) - fft(t.err.di1((1:t.err.L)+t.err.cp1)); %  fft(t.err.do1);
t.err.eq2(t.err.idxSc0) = 1;
t.err.evm(t.err.idxSc0) = 0; 

t.err.do      = t.err.do2 - t.err.do1;
t.err.mse     = sum(abs(t.err.do).^2)/length(t.err.do);

if t.cfg.src == 7,          % OFDM
    x = fft(t.err.do1);
    y = x./t.err.eq1;
    z = fft(t.err.do2);
    z = z./t.err.eq1;
    %scatterplot(y);
    %scatterplot(z);
end

%-----------------------%
%   Plot IQ signal      %
%-----------------------%
if t.iq.plot ==1
    figure(t.iq.figId);
    hold off;
    set(gcf,'Name',t.iq.name);
    subplot(221); 
    hold off;
    stem(real(t.iq.do),'filled');
    hold on;
    line([t.iq.n1 t.iq.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    line([t.iq.n2 t.iq.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    xlim([t.iq.n1,t.iq.n2]);
    xlabel('n'); ylabel('x_I(n)');
    subplot(222);
    plot(freqspace(t.iq.L), mag2db( abs(fftshift(fft(t.iq.do(t.iq.n1:t.iq.n2)))) ));
    xlabel('\Omega/\pi'); ylabel('|X(e^{j\Omega})|   [dB]');
    grid on;
    subplot(223); 
    hold off;
    stem(imag(t.iq.do),'filled');
    hold on;
    line([t.iq.n1 t.iq.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    line([t.iq.n2 t.iq.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    xlim([t.iq.n1,t.iq.n2]);
    xlabel('n'); ylabel('x_Q(n)');
    subplot(224);
    plot(freqspace(t.iq.L), angle(fftshift(fft((t.iq.do(t.iq.n1:t.iq.n2))))));
    xlabel({'\Omega/\pi'}); ylabel('arg\{X(e^{j\Omega})\}');
    grid on;
end

%-----------------------%
% Plot Upsampled signal %
%-----------------------%
if t.rrcFilt.plot == 1
    figure(t.rrcFilt.figId);
    set(gcf,'Name',t.rrcFilt.name);
    subplot(221);
    hold off;
    stem(real(t.rrcFilt.do),'filled');
    hold on;
    line([t.rrcFilt.n1 t.rrcFilt.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    line([t.rrcFilt.n2 t.rrcFilt.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    xlim([t.rrcFilt.n1, t.rrcFilt.n2]);
    xlabel('n'); ylabel('y_I(n)');
    subplot(222);
    plot( freqspace(t.rrcFilt.L), mag2db( abs( fftshift(fft(t.rrcFilt.do(t.rrcFilt.n1:t.rrcFilt.n2))))));
    xlabel('\Omega/\pi'); ylabel('|Y(e^{j\Omega})|  [dB]');
    grid on;
    hold on;
    line([1/t.rrcFilt.sps 1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1/t.rrcFilt.sps -1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    hold off;
    subplot(223);
    hold off;
    stem(imag(t.rrcFilt.do),'filled');
    hold on;
    line([t.rrcFilt.n1 t.rrcFilt.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    line([t.rrcFilt.n2 t.rrcFilt.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    xlim([t.rrcFilt.n1, t.rrcFilt.n2]);
    xlabel('n'); ylabel('y_Q(n)');
    subplot(224);
    plot(freqspace(t.rrcFilt.L), angle(fftshift(fft(t.rrcFilt.do(t.rrcFilt.n1:t.rrcFilt.n2)))));
    xlabel('\Omega/\pi'); ylabel('arg\{Y(e^{j\Omega})\}');
    grid on;
    hold on;
    line([1/t.rrcFilt.sps 1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1/t.rrcFilt.sps -1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    hold off;
end

%-----------------------%
% Plot Filtered signal  %
%-----------------------%
if t.lpFilt.plot == 1
    figure(t.lpFilt.figId);
    set(gcf,'Name',t.lpFilt.name);
    subplot(221);
    hold off;
    stem(real(t.lpFilt.do),'filled');
    hold on;
    line([t.lpFilt.n1 t.lpFilt.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    line([t.lpFilt.n2 t.lpFilt.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    xlim([t.lpFilt.n1, t.lpFilt.n2]);
    xlabel('n'); ylabel('v_I(n)');
    subplot(222);
    plot( freqspace(t.lpFilt.L), mag2db( abs( fftshift(fft(t.lpFilt.do(t.lpFilt.n1:t.lpFilt.n2))))));
    xlabel('\Omega/\pi'); ylabel('|V(e^{j\Omega})|  [dB]');
    grid on;
    hold on;
    line([ 1/t.rrcFilt.sps  1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1/t.rrcFilt.sps -1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    hold off;
    subplot(223);
    hold off;
    stem(imag(t.lpFilt.do),'filled');
    hold on;
    line([t.lpFilt.n1 t.lpFilt.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    line([t.lpFilt.n2 t.lpFilt.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    xlim([t.lpFilt.n1 t.lpFilt.n2]);
    xlabel('n'); ylabel('v_Q(n)');
    subplot(224);
    plot(freqspace(t.lpFilt.L), angle(fftshift(fft(t.lpFilt.do(t.lpFilt.n1:t.lpFilt.n2)))));
    xlabel('\Omega/\pi'); ylabel('arg\{V(e^{j\Omega})\}');
    grid on;
    hold on;
    line([1/t.rrcFilt.sps 1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1/t.rrcFilt.sps -1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    hold off;
end

%-----------------------%
% Plot Theta/Phi signal %
%-----------------------%
if t.angle.plot ==1
    figure(t.angle.figId);
    set(gcf, 'Name', t.angle.name);
    
    subplot(421); 
    hold off;
    stem(t.angle.phi,'filled');
    hold on;
    line([t.angle.n1 t.angle.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    line([t.angle.n2 t.angle.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    xlim([t.angle.n1,t.angle.n2]);
    xlabel('n'); ylabel('\phi(n)');
    subplot(423); 
    hold on; cla;
    
    for i=1:length(t.angle.theta)
        if iscell(t.angle.theta)
            stem(t.angle.theta{i});
        else
            stem(t.angle.theta,'filled');
        end
    end
    
    line([t.angle.n1 t.angle.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    line([t.angle.n2 t.angle.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    xlim([t.angle.n1,t.angle.n2]);
    xlabel('n'); ylabel('\theta(n) ');
    
    subplot(422);
    hold on; cla;
    plot(freqspace(t.angle.L), mag2db(abs(fftshift(fft(t.angle.phi(t.angle.n1:t.angle.n2))))));
    xlabel('\Omega/\pi'); ylabel('|\Phi(e^{j\Omega})|  [dB]');
    grid on;
    line([1/t.rrcFilt.sps 1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1/t.rrcFilt.sps -1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    
    subplot(424); 
    hold on; cla;
    for i=1:length(t.angle.theta)
        if iscell(t.angle.theta)
            plot(freqspace(t.angle.L), mag2db(abs(fftshift(fft(t.angle.theta{i}(t.angle.n1:t.angle.n2))))));
        else
            plot(freqspace(t.angle.L), mag2db(abs(fftshift(fft(t.angle.theta(t.angle.n1:t.angle.n2))))));
        end
    end
    
    xlabel('\Omega/\pi'); ylabel('|\Theta(e^{j\Omega})|  [dB]');
    grid on;
    line([1/t.rrcFilt.sps 1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1/t.rrcFilt.sps -1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    
    subplot(425); 
    hold on; cla;
    if t.cfg.mode ~= 3
        stem(t.angle.alpha,'filled');
    else
        stem(t.angle.theta{1});
        stem(t.angle.theta{2});
    end
    line([t.angle.n1 t.angle.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    line([t.angle.n2 t.angle.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    xlim([t.angle.n1,t.angle.n2]);
    xlabel('n'); ylabel('\phi(n)+\theta(n)');
    
    subplot(427); 
    hold on; cla;
    if t.cfg.mode ~= 3
        stem(t.angle.beta, 'filled');
    else
        stem(-t.angle.theta{1});
        stem(-t.angle.theta{2});
    end
     
    line([t.angle.n1 t.angle.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    line([t.angle.n2 t.angle.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    xlim([t.angle.n1,t.angle.n2]);
    xlabel('n'); ylabel('\phi(n)-\theta(n) ');
    
    subplot(426);
    hold on; cla;
    if t.cfg.mode ~= 3
        plot(freqspace(t.angle.L), mag2db(abs(fftshift(fft(t.angle.alpha(t.angle.n1:t.angle.n2))))));
    else
        plot(freqspace(t.angle.L), mag2db(abs(fftshift(fft(t.angle.theta{1}(t.angle.n1:t.angle.n2) )))));
        plot(freqspace(t.angle.L), mag2db(abs(fftshift(fft(t.angle.theta{2}(t.angle.n1:t.angle.n2) )))));
    end
    xlabel('\Omega/\pi'); ylabel('|(\Phi+\Theta)(e^{j\Omega})|  [dB]');
    grid on;
    line([1/t.rrcFilt.sps 1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1/t.rrcFilt.sps -1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    
    subplot(428); 
    grid on; cla;
    if t.cfg.mode ~= 3
        plot(freqspace(t.angle.L), mag2db(abs(fftshift(fft(t.angle.beta(t.angle.n1:t.angle.n2))))));
    else
        
    end
    xlabel('\Omega/\pi'); ylabel('|(\Phi-\Theta)(e^{j\Omega})|  [dB]');
    hold on;
    line([1/t.rrcFilt.sps 1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1/t.rrcFilt.sps -1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
end

%-----------------------%
% Plot Phase Modulation %
%-----------------------%
if t.pm.plot == 1
    figure(t.pm.figId);
    set(gcf, 'Name', t.pm.name);
    subplot(421);
    hold off;
    stem(real(t.pm.do{1}),'filled');
    hold on;
    line([t.pm.n1 t.pm.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    line([t.pm.n2 t.pm.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    xlim([t.pm.n1, t.pm.n2]);
    xlabel('n'); ylabel('s1_I(n)=sin(\phi(n) + \theta(n))');
    subplot(423);
    hold off;
    stem(imag(t.pm.do{1}),'filled');
    hold on;
    line([t.pm.n1 t.pm.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    line([t.pm.n2 t.pm.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    xlim([t.pm.n1, t.pm.n2]);
    xlabel('n'); ylabel('s1_Q(n)=-cos(\phi(n) + \theta(n))');
    subplot(425);
    hold off;
    stem(real(t.pm.do{2}),'filled');
    hold on;
    line([t.pm.n1 t.pm.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    line([t.pm.n2 t.pm.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    xlim([t.pm.n1, t.pm.n2]);
    xlabel('n'); ylabel('s2_I(n)=sin(\phi(n) - \theta(n))');
    subplot(427);
    hold off;
    stem(imag(t.pm.do{2}),'filled');
    hold on;
    line([t.pm.n1 t.pm.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    line([t.pm.n2 t.pm.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    xlim([t.pm.n1, t.pm.n2]);
    xlabel('n'); ylabel('s2_Q(n)=-cos(\phi(n) - \theta(n))');
    
    subplot(422);
    plot(freqspace(t.pm.L), mag2db((abs(fftshift(fft(t.pm.do{1}(t.pm.n1:t.pm.n2)))))));
    xlabel('\Omega/\pi'); ylabel('|S_I(e^{j\Omega}|  [dB]');
    grid on;
    hold on;
    line([1/t.rrcFilt.sps 1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1/t.rrcFilt.sps -1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    hold off;
    subplot(424);
    plot(freqspace(t.pm.L), angle(fftshift(fft(t.pm.do{1}(t.pm.n1:t.pm.n2)))));
    xlabel('\Omega/\pi'); ylabel('arg\{S_I(e^{j\Omega}\}');
    grid on;
    hold on;
    line([1/t.rrcFilt.sps 1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1/t.rrcFilt.sps -1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    hold off;
    
    subplot(426);
    plot(freqspace(t.pm.L), mag2db((abs(fftshift(fft(t.pm.do{2}(t.pm.n1:t.pm.n2)))))));
    xlabel('\Omega/\pi'); ylabel('|S_Q(e^{j\Omega}|  [dB]');
    grid on;
    hold on;
    line([1/t.rrcFilt.sps 1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1/t.rrcFilt.sps -1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    hold off;
    subplot(428);
    plot(freqspace(t.pm.L), angle(fftshift(fft(t.pm.do{2}(t.pm.n1:t.pm.n2)))));
    xlabel('\Omega/\pi'); ylabel('arg\{S_Q(e^{j\Omega}\}');
    grid on;
    hold on;
    line([1/t.rrcFilt.sps 1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1/t.rrcFilt.sps -1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    hold off;
end

if t.aFilt.plot == 1
    figure(t.aFilt.figId);
    set(gcf, 'Name', t.aFilt.name);
    
    opt = {'ai','bi'; 'a','b'};       % input coef field names

    for m=1:2,
        for i=1:2,
        
            if m == 1
                ref = 1;
            else
                ref = 0;
            end
        
            a = t.aFilt.(opt{m,1});
            b = t.aFilt.(opt{m,2});
        
            subplot(4,2, 4*(m-1) + 2*i-1);
            W = linspace( -(t.aFilt.wc+t.aFilt.bw), t.aFilt.wc+t.aFilt.bw ,500);
            %H = freqs(t.aFilt.bi{i},t.aFilt.ai{i}, W);
            H = freqs(b{i},a{i}, W);
            plot(W/t.aFilt.wc, mag2db(abs(H)));
            grid on;
            xlabel('\omega/\omega_c');
            ylabel('Magnitude (dB)');
            hold on;
            line([ref+t.aFilt.bw/2/t.aFilt.wc ref+t.aFilt.bw/2/t.aFilt.wc],get(gca,'YLim'),'Color',[1 0 0]);
            line([ref-t.aFilt.bw/2/t.aFilt.wc ref-t.aFilt.bw/2/t.aFilt.wc],get(gca,'YLim'),'Color',[1 0 0]);
            hold off;
    
            subplot(4,2, 4*(m-1) + 2*i);
            plot(W/t.aFilt.wc,phase(H));
            xlabel('\omega/\omega_c');
            ylabel('Phase');
            grid on;
            hold on;
            line([ref+t.aFilt.bw/2/t.aFilt.wc ref+t.aFilt.bw/2/t.aFilt.wc],get(gca,'YLim'),'Color',[1 0 0]);
            line([ref-t.aFilt.bw/2/t.aFilt.wc ref-t.aFilt.bw/2/t.aFilt.wc],get(gca,'YLim'),'Color',[1 0 0]);
            hold off;
        end
    end
end

%-----------------------%
% Plot Baseband Filter  %
%-----------------------%
if t.bbFilt.plot ==1
    figure(t.bbFilt.figId);
    set(gcf, 'Name', t.bbFilt.name);
    
    
    for i=1:2,      
        a = t.bbFilt.a{i};
        b = t.bbFilt.b{i};
        
        W = linspace( -pi, pi ,500);
        [H] = freqz(b, a, W);
        
        subplot(6,2, 2*i-1 );
        plot(W/pi,mag2db(abs(H)));
        grid on;
        xlabel('\Omega/\pi');
        ylabel('Magnitude [dB]');
        hold on;
        line([ 1/t.rrcFilt.sps  1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
        line([-1/t.rrcFilt.sps -1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
        hold off;
        
        subplot(6,2, 2*i);
        plot(W/pi,angle(H));
        grid on;
        xlabel('\Omega/\pi');
        ylabel('Phase');
        hold on;
        line([ 1/t.rrcFilt.sps  1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
        line([-1/t.rrcFilt.sps -1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
        hold off;
    end
    
    % Output signal
    for i=1:2,
        do=t.bbFilt.do{i};
    
        subplot(6,2, 5+4*(i-1));
        hold off;
        stem(real(do),'filled');
        hold on;
        line([t.bbFilt.n1 t.bbFilt.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
        line([t.bbFilt.n2 t.bbFilt.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
        xlim([t.bbFilt.n1, t.bbFilt.n2]);
        xlabel('n'); ylabel('Re(n)');
        
        subplot(6,2,5+4*(i-1)+2);
        hold off;
        stem(imag(do),'filled');
        hold on;
        line([t.bbFilt.n1 t.bbFilt.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
        line([t.bbFilt.n2 t.bbFilt.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
        xlim([t.bbFilt.n1 t.bbFilt.n2]);
        xlabel('n'); ylabel('Im(n)');
        
        
        subplot(6,2, 6+4*(i-1));
        plot( freqspace(t.bbFilt.L), mag2db( abs( fftshift(fft(do(t.bbFilt.n1:t.bbFilt.n2))))));
        xlabel('\Omega/\pi'); ylabel('Magnitude  [dB]');
        grid on;
        hold on;
        line([ 1/t.rrcFilt.sps  1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
        line([-1/t.rrcFilt.sps -1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
        hold off;
        
        
        subplot(6,2, 6+4*(i-1)+2);
        plot(freqspace(t.bbFilt.L), angle(fftshift(fft(do(t.bbFilt.n1:t.bbFilt.n2)))));
        xlabel('\Omega/\pi'); ylabel('Phase');
        grid on;
        hold on;
        line([1/t.rrcFilt.sps 1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
        line([-1/t.rrcFilt.sps -1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
        hold off;
    end
    
end

%-----------------------%
% Plot Combiner Signal  %
%-----------------------%
if t.comb.plot == 1
    figure(t.comb.figId);
    set(gcf,'Name',t.comb.name);
    
    subplot(221);
    hold off;
    stem(real(t.comb.do),'filled');
    hold on;
    line([t.comb.n1 t.comb.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    line([t.comb.n2 t.comb.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    xlim([t.comb.n1, t.comb.n2]);
    xlabel('n'); ylabel('c_I(n)');
    subplot(222);
    plot( freqspace(t.comb.L), mag2db( abs( fftshift(fft( t.comb.do(t.comb.n1:t.comb.n2))))));
    xlabel('\Omega/\pi'); ylabel('|C(e^{j\Omega})|  [dB]');
    grid on;
    hold on;
    line([1/t.rrcFilt.sps 1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1/t.rrcFilt.sps -1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    hold off;
    
    subplot(223);
    hold off;
    stem(imag(t.comb.do),'filled');
    hold on;
    line([t.comb.n1 t.comb.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    line([t.comb.n2 t.comb.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    xlim([t.comb.n1, t.comb.n2]);
    xlabel('n'); ylabel('c_Q(n)');
    subplot(224);
    plot(freqspace(t.comb.L), angle(fftshift(fft(t.comb.do(t.comb.n1:t.comb.n2)))));
    xlabel('\Omega/\pi'); ylabel('arg\{C(e^{j\Omega})\}');
    grid on;
    hold on;
    line([1/t.rrcFilt.sps 1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1/t.rrcFilt.sps -1/t.rrcFilt.sps],get(gca,'YLim'),'Color',[1 0 0]);
    hold off;
end

%--------------------------%
% Plot Downsampled Signal  %
%--------------------------%
if t.irrcFilt.plot == 1
    figure(t.irrcFilt.figId);
    set(gcf,'Name',t.irrcFilt.name);
    
    subplot(221);
    hold off;
    stem(real(t.irrcFilt.do),'filled');
    hold on;
    line([t.irrcFilt.n1 t.irrcFilt.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    line([t.irrcFilt.n2 t.irrcFilt.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    xlim([t.irrcFilt.n1, t.irrcFilt.n2]);
    xlabel('n'); ylabel('y_I(n)');
    subplot(222);
    plot( freqspace(t.irrcFilt.L), mag2db( abs( fftshift(fft(t.irrcFilt.do(t.irrcFilt.n1:t.irrcFilt.n2))))));
    xlabel('\Omega/\pi'); ylabel('|Y(e^{j\Omega})|  [dB]');
    grid on;
    hold on;
    line([1 1],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1 -1],get(gca,'YLim'),'Color',[1 0 0]);
    hold off;
    subplot(223);
    hold off;
    stem(imag(t.irrcFilt.do),'filled');
    hold on;
    line([t.irrcFilt.n1 t.irrcFilt.n1],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    line([t.irrcFilt.n2 t.irrcFilt.n2],get(gca,'YLim'),'Color',[1 0 0],'LineStyle','--');
    xlim([t.irrcFilt.n1, t.irrcFilt.n2]);
    xlabel('n'); ylabel('y_Q(n)');
    subplot(224);
    plot(freqspace(t.irrcFilt.L), angle(fftshift(fft(t.irrcFilt.do(t.irrcFilt.n1:t.irrcFilt.n2)))));
    xlabel('\Omega/\pi'); ylabel('arg\{Y(e^{j\Omega})\}');
    grid on;
    hold on;
    line([1 1],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1 -1],get(gca,'YLim'),'Color',[1 0 0]);
    hold off;
end

%--------------------%
% Plot Error Signal  %
%--------------------%
if t.err.plot == 1
    figure(t.err.figId);
    set(gcf,'Name',t.err.name);
    
    subplot(421);
    hold off;
    stem(real(  t.err.do1) , 'filled');
    hold on;
    stem(real( t.err.do2),'r');
    xlabel('n'); ylabel('Re\{In\}, Re\{Out\} (red)');
    subplot(422);
    hold off;
    stem(imag(  t.err.do1), 'filled');
    hold on;
    stem(imag( t.err.do2 ),'r');
    xlabel('n'); ylabel('Im\{In\}, Im\{Out\} (red)');
    
    
    subplot(423);
    stem(abs(t.err.do),'filled');
    xlabel('n'); ylabel(t.err.ylabel);
    title(['MSE = ', num2str(t.err.mse)]);
    subplot(424);
    stem(angle(t.err.do),'filled');
    xlabel('n'); ylabel('arg\{e(n)\}');
    
    subplot(425);
    plot( freqspace(length(t.err.do)), mag2db( abs( fftshift(fft(t.err.do)))));
    xlabel('\Omega/\pi'); ylabel('|E(e^{j\Omega})|  [dB]');
    grid on;
    hold on;
    line([ 1/t.err.sps  1/t.err.sps],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1/t.err.sps -1/t.err.sps],get(gca,'YLim'),'Color',[1 0 0]);
    hold off;
    subplot(426);
    plot(freqspace(length(t.err.do)), angle(fftshift(fft(t.err.do))));
    xlabel('\Omega/\pi'); ylabel('arg\{E(e^{j\Omega})\}');
    grid on;
    hold on;
    line([ 1/t.err.sps  1/t.err.sps],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1/t.err.sps -1/t.err.sps],get(gca,'YLim'),'Color',[1 0 0]);
    hold off;
    
    subplot(427);
    plot( freqspace(length(t.err.do)),  abs( fftshift(t.err.eq2)));
    xlabel('\Omega/\pi'); ylabel('|H(e^{j\Omega})|');
    grid on;
    hold on;
    line([ 1/t.err.sps  1/t.err.sps],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1/t.err.sps -1/t.err.sps],get(gca,'YLim'),'Color',[1 0 0]);
    hold off;
    subplot(428);
    plot(freqspace(length(t.err.do)), abs(fftshift(t.err.evm)));
    xlabel('\Omega/\pi'); ylabel('|EVM(e^{j\Omega})| ');
    grid on;
    hold on;
    line([ 1/t.err.sps  1/t.err.sps],get(gca,'YLim'),'Color',[1 0 0]);
    line([-1/t.err.sps -1/t.err.sps],get(gca,'YLim'),'Color',[1 0 0]);
    hold off;
end
%}

end % end function



%--------------------------------------------%
%  subfunctions                              %
%--------------------------------------------%


function diff = delay(x,y),
% Compute delay of x and y vectors

[c, lags] = xcorr(x, y);
[~,I] = max(abs(c));
diff = lags(I);

end

function y = roots2z(x, fs),
% Transforms pole or zero from s domain to z domain using bilinear
% transformation
%
% x  - pole/zero in s-domain
% fs - sampling frequency in Hz
% y  - pole/zero in z domain
y = (2*fs + x)/(2*fs - x);
end

function y = freqspace(L),
% Compute linear spaced L-length vector of normalized freq. bins in interval <-1,1)

    y = (-L/2:L/2-1)/L*2;

end
