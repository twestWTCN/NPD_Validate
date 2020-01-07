function [f,t,cl,fyxzw]=sp2_R2a_pc1_tw(x,y,z,samp_rate,seg_size,phrand)
if nargin<6
    phrand = 0;
end
% function [f,t,cl,fyxw]=sp2_R2a_pc1(x,y,z,samp_rate,seg_size)
%
% Two channel average periodogram analysis, which include R2 analysis (ver 2) using MMSE pre-whitening
%  which generates fxx and fyy identically 1 at all frequencies. This version for fyx/z, 1st order partial.
%
% Only pre-processing is mean subtraction.
%
% Output is sp2* format f, t, (extra columns in f and t) with different cl structure
% Inputs
%       x          Vector with random process x
%       y          Vector with random process y
%       z          Vector with random process z, predictor
%       samp_rate  Sampling rate - to get Fourier frequencies
%       seg_size   Segment size to use in analysis
%
% Outputs
%       f, t, cl, fyxzw

% Create zero mean sequence
x=x-mean(x);
y=y-mean(y);
z=z-mean(z);

% Segment for FFT - rows: seg_size. Columns: seg_tot
seg_tot=floor(length(x)/seg_size);
samp_tot=seg_tot*seg_size;
x=reshape(x(1:seg_tot*seg_size),seg_size,seg_tot);
y=reshape(y(1:seg_tot*seg_size),seg_size,seg_tot);
z=reshape(z(1:seg_tot*seg_size),seg_size,seg_tot);

% Zero mean for each segment
for ind=1:seg_tot
	x(:,ind)=x(:,ind)-mean(x(:,ind));
	y(:,ind)=y(:,ind)-mean(y(:,ind));
	z(:,ind)=z(:,ind)-mean(z(:,ind));
end

% FFT of 3 processes
dx=fft(x);
dy=fft(y);
dz=fft(z);

if phrand
    phaseRandSym(dx)
    
end

% Need fxz, fyz, fzz to generate conditional dFT
Izz=1/(2*pi*seg_size)*abs(dz.*dz);
Ixz=1/(2*pi*seg_size)*dx.*conj(dz);
Iyz=1/(2*pi*seg_size)*dy.*conj(dz);

fzz=mean(Izz,2);
fxz=mean(Ixz,2);
fyz=mean(Iyz,2);

% Gain factors to generate conditional dFT.	26/2/16 DC forced to zero
%gxz=fxz./fzz;
%gyz=fyz./fzz;
gxz=[0; fxz(2:seg_size)./fzz(2:seg_size)];
gyz=[0; fyz(2:seg_size)./fzz(2:seg_size)];

% Generate conditional dFTs: dx|z, dy|z
dxz=dx;
dyz=dy;
for ind=1:seg_tot
	dxz(:,ind)=dxz(:,ind)-gxz.*dz(:,ind);
	dyz(:,ind)=dyz(:,ind)-gyz.*dz(:,ind);
end

% Partial spectra and coherence from conditional dFT
Ixxz=1/(2*pi*seg_size)*abs(dxz.*dxz);
Iyyz=1/(2*pi*seg_size)*abs(dyz.*dyz);
Iyxz=1/(2*pi*seg_size)*dyz.*conj(dxz);

fxxz=mean(Ixxz,2);
fyyz=mean(Iyyz,2);
fyxz=mean(Iyxz,2);
%chyxz=abs(fyxz.*fyxz)./(fxxz.*fyyz);  % 11/12/13 DC forced to zero
chyxz=[0; abs(fyxz(2:seg_size).*fyxz(2:seg_size))./(fxxz(2:seg_size).*fyyz(2:seg_size))];

fxx=fxxz;
fyy=fyyz;
fyx=fyxz;
chyx=chyxz;

%Pre-whitening stage - generate MMSE filters	11/12/13 weight for DC forced to zero
%wx=1./sqrt(fxxz);
%wy=1./sqrt(fyyz);
wx=[0; 1./sqrt(fxxz(2:seg_size))];
wy=[0; 1./sqrt(fyyz(2:seg_size))];

% Apply MMSE pre-whitening filters
dxzw=dxz;
dyzw=dyz;
for ind=1:seg_tot
	dxzw(:,ind)=dxzw(:,ind).*wx;
	dyzw(:,ind)=dyzw(:,ind).*wy;
end

% Re-calculate using pre-whitened processes
% Using bivariate variable names 
Ixxw=1/(2*pi*seg_size)*abs(dxzw.*dxzw);
Iyyw=1/(2*pi*seg_size)*abs(dyzw.*dyzw);
Iyxw=1/(2*pi*seg_size)*dyzw.*conj(dxzw);
fxxw=mean(Ixxw,2);      % Should be 1 all freqs
fyyw=mean(Iyyw,2);      %        "
fyxw=mean(Iyxw,2);      % Complex cross spectrum
chyxw=abs(fyxw.*fyxw);  % Coherence now from fyxw only
rhoyx=real(ifft(fyxw)); % rhoyx, NB no factors, 1/T in ifft()

% R2 values
R2_ch=(1/seg_size)*sum(chyx);   % Integral from -pi to +pi
R2_chw=(1/seg_size)*sum(chyxw); % Integral from -pi to +pi

% R2 in time domain
R2_rho=sum(rhoyx.^2);
% Decompose R2 by lag
R2_rho_0=rhoyx(1).^2;                          % zero lag, bin 1
R2_rho_p=sum(rhoyx(2:seg_size/2).^2);          % positive lags, bins 2:(T/2)
R2_rho_n=sum(rhoyx(seg_size/2+1:seg_size).^2); % negative lags, bins (T/2)+1:T

% Decompose rho into three components by lag
rhoyx_3=zeros(length(rhoyx),3);
rhoyx_3(1,1)=rhoyx(1);                                         % zero lag
rhoyx_3(2:seg_size/2,2)=rhoyx(2:seg_size/2);                   % positive lags
rhoyx_3(seg_size/2+1:seg_size,3)=rhoyx(seg_size/2+1:seg_size); % negative lags

% Switch back to frequency domain
f_prime=fft(rhoyx_3);      % These are f' values:     f'yx;0, f'yx;+, f'yx;-
f_prime2=abs(f_prime).^2;  % These are |f'|^2 values: |f'yx;0|^2, |f'yx;+|^2, |f'yx;-|^2

% Scale to get relative contributions from 0,+,- components
% These will have range 0<= ch_3 <=1, these are multipliers in (current) eqns (20)-(22)
f_prime2_fac=sum(f_prime2,2);
R2_weight=f_prime2;
R2_weight(:,1)=R2_weight(:,1)./f_prime2_fac;
R2_weight(:,2)=R2_weight(:,2)./f_prime2_fac;
R2_weight(:,3)=R2_weight(:,3)./f_prime2_fac;

% R^2 decomposition using integrated |f'|^2 - should agree with decomposition by lag
R2_fprime_0=(1/seg_size)*sum(f_prime2(:,1));
R2_fprime_p=(1/seg_size)*sum(f_prime2(:,2));
R2_fprime_n=(1/seg_size)*sum(f_prime2(:,3));
% R2 decomposition using |R'|^2 values: |R'yx;0|^2, |R'yx;+|^2, |R'yx;-|^2
R2_Rprime_0=(1/seg_size)*sum(R2_weight(:,1).*chyx);
R2_Rprime_p=(1/seg_size)*sum(R2_weight(:,2).*chyx);
R2_Rprime_n=(1/seg_size)*sum(R2_weight(:,3).*chyx);

%-----------------------------------------------------------------------
% Output matrix f
f=zeros(seg_size/2,15);
f_index=(2:seg_size/2+1)';    % Indexing for output, DC not output.
deltaf=samp_rate/seg_size;
f(:,1)=(f_index-1)*deltaf;    % Column 1 - frequencies in Hz.
f(:,2)=log10(fxx(f_index));   % Column 2 - Log spectrum ch 1.
f(:,3)=log10(fyy(f_index));   % Column 3 - Log spectrum ch 2.
f(:,4)=chyx(f_index);         % Column 4 - Coherence, PBMB (5.5).
f(:,5)=angle(fyx(f_index));   % Column 5 - Phase, PBMB (5.7).
% Append to f - spectra NOT log10, should be 1 all freqs, coherence from fyxw only.
f(:,6)=fxxw(f_index);
f(:,7)=fyyw(f_index);
f(:,8)=chyxw(f_index);
f(:,9)=angle(fyxw(f_index));
f(:,10)=R2_weight(f_index,1).*chyx(f_index);   % Component from zero lag (will be const)
f(:,11)=R2_weight(f_index,2).*chyx(f_index);   % Directional component positive lag, R2+
f(:,12)=R2_weight(f_index,3).*chyx(f_index);   % Directional component negative lag, R2-
f(:,13)=f_prime2(f_index,1);   % f'yx;0
f(:,14)=f_prime2(f_index,2);   % f'yx;+
f(:,15)=f_prime2(f_index,3);   % f'yx;-

% Time domain output: Col 1 - lag  range (-T/2)*dt to (T/2-1)*dt, Col 2 - cumulant, Col 4 - rho
deltat=1000.0/samp_rate;    % dt in msec.
t(:,1)=((1:seg_size)'-seg_size/2-1)*deltat;
t([seg_size/2+1:seg_size,1:seg_size/2],2)=2*pi*real(ifft(fyx));
t([seg_size/2+1:seg_size,1:seg_size/2],3)=rhoyx;  % Col 3 is rhoyx, same indexing as q in col 2
% Estimate variance of cumulant density estimate.
var_fac=4*pi*pi/(seg_size*samp_tot);                           % Factor (2pi/S)(2pi/R).
q_var=var_fac*2*sum(fxx(1:seg_size/2+1).*fyy(1:seg_size/2+1)); % PBMB (6.10).
rho_var=1/samp_tot;

cl.type=0;                   % Analysis type.
cl.seg_size=seg_size;        % S.
cl.seg_tot=seg_tot;          % L.
cl.seg_tot_var=seg_tot;      % Effective no of segments (L').
cl.samp_tot=seg_tot*seg_size;% R.
cl.samp_rate=samp_rate;      % Sampling rate.
cl.dt=deltat;                % Delta t.
cl.df=deltaf;                % Delta f.
cl.f_c95=0.8512*sqrt(1/cl.seg_tot_var);  % 95% Confidence limit for spectral estimates, PBMB (6.2).
% N.B. Confidence interval for log plot of spectra is TWICE this value.
cl.ch_c95=1-0.05^(1/(cl.seg_tot_var-1)); % 95% Confidence limit for coherence, PBMB (6.6).
cl.q_c95=1.96*sqrt(q_var);               % 95% Confidence limits for cumulant, PBMB (6.11).
cl.rho_c95=1.96*sqrt(rho_var);
cl.col_R20=10;               % Cols for plotting
cl.col_R2p=11;
cl.col_R2n=12;
cl.col_rho=3;
cl.R2_ch=R2_ch;
cl.R2_chw=R2_chw;
cl.R2_rho=R2_rho;
cl.R2_rho_0=R2_rho_0;
cl.R2_rho_n=R2_rho_n;
cl.R2_rho_p=R2_rho_p;
cl.R2_fprime_0=R2_fprime_0;
cl.R2_fprime_n=R2_fprime_n;
cl.R2_fprime_p=R2_fprime_p;
cl.R2_Rprime_0=R2_Rprime_0;
cl.R2_Rprime_n=R2_Rprime_n;
cl.R2_Rprime_p=R2_Rprime_p;
cl.R2_ch_dc=chyx(1);
cl.R2_fprime2_dc=[f_prime2(1,1) f_prime2(1,2) f_prime2(1,3)];
cl.R2_Rprime2_dc=[R2_weight(1,1).*chyx(1) R2_weight(1,2).*chyx(1) R2_weight(1,3).*chyx(1)];
cl.N2=0;        % N2, No of events analysed from sp1.
cl.N1=0;        % N1, No of events analysed from sp2.
cl.P2=0;        % P2, mean intensity of events from sp1.
cl.P1=0;        % P1, mean intensity of events from sp2.
cl.opt_str='';  % Options string.
cl.what='';     % Field for plot label.

% Display messages
%disp(['Segments: ',num2str(seg_tot)])

disp(['R2 (R2_ch, R2_chw, R2_rho): ',num2str([R2_ch, R2_chw, R2_rho])])
disp(['  R2_rho values (0, n, p)): ',num2str([R2_rho_0, R2_rho_n, R2_rho_p]),'.  Sum: ',num2str([R2_rho_0+R2_rho_n+R2_rho_p])])
disp(['  R2_f''  values (0, n, p)): ',num2str([R2_fprime_0, R2_fprime_n, R2_fprime_p]),'.  Sum: ',num2str([R2_fprime_0+R2_fprime_n+R2_fprime_p])])
disp(['  R2_R''  values (0, n, p)): ',num2str([R2_Rprime_0, R2_Rprime_n, R2_Rprime_p]),'.  Sum: ',num2str([R2_Rprime_0+R2_Rprime_n+R2_Rprime_p])])

%disp([' ch diff: ',num2str([min(f(:,4)-f(:,8)) max(f(:,4)-f(:,8))])])
%disp([' ph diff: ',num2str([min(f(:,5)-f(:,9)) max(f(:,5)-f(:,9))])])
