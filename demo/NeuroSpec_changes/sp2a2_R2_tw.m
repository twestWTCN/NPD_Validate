function [f,t,cl]=sp2a2_R2_tw(x,y,samp_rate,seg_pwr)
% function [f,t,cl]=sp2a2_R2(x,y,samp_rate,seg_pwr)
%
% Two channel average periodogram analysis, with R2 directionality analysis.
%
% Copyright 2015, David M. Halliday.
% This file is part of NeuroSpec.
%
%    NeuroSpec is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    NeuroSpec is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with NeuroSpec; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%    NeuroSpec is available at:  http://www.neurospec.org/
%    Contact:  contact@neurospec.org
%
%
% Input arguments
%  x          Channel 1  (input) time series vector
%  y          Channel 2 (output) time series vector
%  samp_rate  Sampling rate (samples/sec)
%  seg_pwr    Segment length - specified as power of 2
%
%
% Output arguments
%  f    column matrix with frequency domain parameter estimates
%  t    column matrix with      time domain parameter estimates
%  cl   single structure with scalar values related to analysis
%
% Output parameters
%  f column 1       frequency in Hz
%  f column 2       Log input/x  spectrum
%  f column 3       Log output/y spectrum
%  f column 4       Coherence
%  f column 5       Phase
%  f column 6       MMSE whitened spectra for process x, JIN (2.4)
%  f column 7       MMSE whitened spectra for process y, JIN (2.4)
%  f column 8       Coherence between whitened processes, JIN (2.5)
%  f column 9       Phase between whitened processes
%  f column 10      Zero lag coherence component, JIN (2.19)
%  f column 11      Forward  coherence component, JIN (2.20)
%  f column 12      Reverse  coherence component, JIN (2.18)
%
%  t column 1       Lag in ms.
%  t column 2       Cumulant/covariance density.
%  t column 3       rhoyx, JIN (2.8)
%
%  cl.type          Analysis type 0, for compatability with NeuroSpec2.0
%  cl.seg_size      Segment length, T
%  cl.seg_tot       Number of segments, L
%  cl.seg_tot_var   Effective no of segments, used to calculate confidence limits
%  cl.samp_tot      Number of samples analysed, R=LT
%  cl.samp_rate     Sampling rate of data (samples/sec)
%  cl.dt            Time domain bin width (ms)
%  cl.df            Frequency domain bin width (Hz)
%  cl.f_c95         95% confidence limit for Log spectral estimates
%  cl.ch_c95        95% confidence limit for coherence estimate
%  cl.q_c95         95% confidence limit for cumulant density estimate
%  cl.rho_c95       95% confidence limit for rho estimate
%  cl.col_R20       column in f matrix containing zero lag coherence component
%  cl.col_R2p       column in f matrix containing forward (positive lag) coherence component
%  cl.col_R2n       column in f matrix containing reverse (negative lag) coherence component
%  cl.col_rho       column in t matrix containing rho estimate
%  cl.R2            R2 value, R^2_yx
%  cl.R2_0          Component of R2 at zero lag, R^2_yx;0
%  cl.R2_p          Component of R2 in forward direction, positive lag component, R^2_yx;+
%  cl.R2_n          Component of R2 in reverse direction, negative lag component, R^2_yx;-
%  cl.R2_ch         R2 estimated from coherence between original processes
%  cl.R2_chw        R2 estimated from coherence between MMSE whitened processes
%  cl.R2_fprime_0   Component of R2 at zero lag estimated from f'
%  cl.R2_fprime_p   Component of R2 at positive lag estimated from f'
%  cl.R2_fprime_n   Component of R2 at negative lag estimated from f'
%  cl.N1, cl.N2     Spike counts for point process data
%  cl.P1, cl.P2     Spike rates for point process data
%  cl.opt_str       Copy of options string, empty here
%  cl.what          Text label, used in plotting routines
%
% References:
% 1. Halliday D.M., Rosenberg J.R., Amjad A.M., Breeze P., Conway B.A. & Farmer S.F.
%     Progress in Biophysics and molecular Biology, 64, 237-278, 1995.
% 2. Halliday DM (2015) Nonparametric directionality measures for time series and point process data,
%     Journal of Integrative Neuroscience, 14(2), In Press. DOI: 10.1142/S0219635215300127
%
% These references referred to in comments as PBMB and JIN, respectively
%
% function [f,t,cl]=sp2a2_R2(x,y,samp_rate,seg_pwr)

% Check number of input arguments
if (nargin<4)
  error(' Not enough input arguments');
end
% Check number of output arguments
if (nargout<3)
  error(' Not enough output arguments');
end

% Check for single column data
[nrow,ncol]=size(x);
if (ncol~=1)
  error(' Input NOT single column: x')
end
[nrow,ncol]=size(y);
if (ncol~=1)
  error(' Input NOT single column: y')
end

pts_tot=length(x);           % Determine size of data vector
if (length(y)~=pts_tot)      % Check that input vectors are equal length
  error (' Unequal length data arrays');
end

if (max(size(samp_rate)) ~= 1)
  error(' Non scalar value for: samp_rate');
end
if (max(size(seg_pwr)) ~= 1)
  error(' Non scalar value for: seg_pwr');
end

% Segment for FFT - Rows: seg_size, T. Columns: seg_tot, L
seg_size=fix(2^seg_pwr);              % Segment length, T
if rem(seg_size,2)
    seg_size = seg_size+1;
end
seg_tot=fix(length(x)/seg_size); % No of segments, L
samp_tot=seg_tot*seg_size;       % Record length,  R=LT
x=reshape(x(1:samp_tot),seg_size,seg_tot);  % T rows, L columns
y=reshape(y(1:samp_tot),seg_size,seg_tot);  % T rows, L columns

% Create zero mean sequence for each segment.
for ind=1:seg_tot
  x(:,ind)=x(:,ind)-mean(x(:,ind));
  y(:,ind)=y(:,ind)-mean(y(:,ind));
end

% FFT & average periodogram
dx=fft(x);   % dFT across columns, PBMB (4.1)
dy=fft(y);

psd_fac=1/(2*pi*samp_tot);
fxx=psd_fac*sum(abs(dx.*dx),2);  % PBMB (5.2)
fyy=psd_fac*sum(abs(dy.*dy),2);  % PBMB (5.2)
fyx=psd_fac*sum(dy.*conj(dx),2); % PBMB (5.2)

% Coherence PBMB (5.4), with zero frequency set to zero
chyx=[0; abs(fyx(2:seg_size).*fyx(2:seg_size))./(fxx(2:seg_size).*fyy(2:seg_size))];

%Pre-whitening stage - generate MMSE filters, weight for zero frequency set to zero
wx=[0; 1./sqrt(fxx(2:seg_size))];  % JIN (2.23)
wy=[0; 1./sqrt(fyy(2:seg_size))];  % JIN (2.24)

% Apply MMSE pre-whitening filters
dxw=dx;
dyw=dy;
for ind=1:seg_tot
  dxw(:,ind)=dxw(:,ind).*wx;  % JIN (2.25)
  dyw(:,ind)=dyw(:,ind).*wy;  % JIN (2.26)
end

% Re-calculate spectra using pre-whitened processes
fxxw=psd_fac*sum(abs(dxw.*dxw),2);   % Spectra of MMSE whitened process x. This is 1 at all freqs
fyyw=psd_fac*sum(abs(dyw.*dyw),2);   % Spectra of MMSE whitened process y. This is 1 at all freqs
fyxw=psd_fac*sum(dyw.*conj(dxw),2);  % Complex cross spectrum between whitened processes
chyxw=abs(fyxw.*fyxw);               % Coherence from fyxw only, JIN (2.27) 
rhoyx=real(ifft(fyxw));              % rhoyx, no factors as 1/T in ifft(), estimate of JIN (2.8) 

% R2 in frequency domain
R2_ch=(1/seg_size)*sum(chyx);   % Integral from -pi to +pi, estimate for JIN (2.2)
R2_chw=(1/seg_size)*sum(chyxw); % Integral from -pi to +pi, JIN (2.28)

% Estimate R2 in time domain and separate components:
% zero lag R^2_yx;0, positive lag R^2_yx;+ and negative lag R^2_yx;-
R2_rho=sum(rhoyx.^2);                          % JIN (2.29)
R2_rho_0=rhoyx(1).^2;                          % JIN (2.31) zero lag, bin 1
R2_rho_p=sum(rhoyx(2:seg_size/2).^2);          % JIN (2.32) positive lags, bins 2:(T/2)
R2_rho_n=sum(rhoyx(seg_size/2+1:seg_size).^2); % JIN (2.30) negative lags, bins (T/2)+1:T

% Decompose rho into three components by lag, to estimate f' functions
rhoyx_3=zeros(length(rhoyx),3);
rhoyx_3(1,1)=rhoyx(1);                                         % zero lag
rhoyx_3(2:seg_size/2,2)=rhoyx(2:seg_size/2);                   % positive lags
rhoyx_3(seg_size/2+1:seg_size,3)=rhoyx(seg_size/2+1:seg_size); % negative lags

% Switch back to frequency domain, these are f' estimates
f_prime=fft(rhoyx_3);      % JIN (2.34, 2.35, 2.33):  f'yx;0, f'yx;+, f'yx;-
f_prime2=abs(f_prime).^2;  % These are |f'|^2 values: |f'yx;0|^2, |f'yx;+|^2, |f'yx;-|^2

% Scale to get relative contributions from 0,+,- components.
f_prime2_fac=sum(f_prime2,2);
R2_weight=f_prime2;
R2_weight(:,1)=R2_weight(:,1)./f_prime2_fac;  % Weighting factor in JIN (2.19), range [0 1]
R2_weight(:,2)=R2_weight(:,2)./f_prime2_fac;  % Weighting factor in JIN (2.20), range [0 1]
R2_weight(:,3)=R2_weight(:,3)./f_prime2_fac;  % Weighting factor in JIN (2.18), range [0 1]

% R2 metrics using integrated |f'|^2 - should agree with decomposition by lag
R2_fprime_0=(1/seg_size)*sum(f_prime2(:,1));  % These are terms in JIN (2.15)
R2_fprime_p=(1/seg_size)*sum(f_prime2(:,2));
R2_fprime_n=(1/seg_size)*sum(f_prime2(:,3));

%-----------------------------------------------------------------------
% Construct output matrix f
f=zeros(seg_size/2,12);
f_index=(2:seg_size/2+1)';    % Indexing for output, DC not output
deltaf=samp_rate/seg_size;
f(:,1)=(f_index-1)*deltaf;    % Column 1 - frequencies in Hz
f(:,2)=log10(fxx(f_index));   % Column 2 - Log spectrum ch 1, x
f(:,3)=log10(fyy(f_index));   % Column 3 - Log spectrum ch 2, y
f(:,4)=chyx(f_index);         % Column 4 - Coherence, PBMB (5.5)
f(:,5)=angle(fyx(f_index));   % Column 5 - Phase, PBMB (5.7)
f(:,6)=fxxw(f_index);         % Column 6 - MMSE whitened spectra process x, JIN (2.4)
f(:,7)=fyyw(f_index);         % Column 7 - MMSE whitened spectra process y, JIN (2.4)
f(:,8)=chyxw(f_index);        % Column 8 - Coherence between whitened processes, JIN (2.5)
f(:,9)=angle(fyxw(f_index));  % Column 9 - Phase between whitened processes
f(:,10)=R2_weight(f_index,1).*chyxw(f_index);   % Column 10 - Zero lag coherence component, JIN (2.19)
f(:,11)=R2_weight(f_index,2).*chyxw(f_index);   % Column 11 - Forward  coherence component, JIN (2.20)
f(:,12)=R2_weight(f_index,3).*chyxw(f_index);   % Column 12 - Reverse  coherence component, JIN (2.18)

% Construct time domain matrix t
% Time domain output: Col 1 - lag  range (-T/2)*dt to (T/2-1)*dt, Col 2 - cumulant, Col 3 - rho
deltat=1000.0/samp_rate;                                         % dt in msec.
t(:,1)=((1:seg_size)'-seg_size/2-1)*deltat;                      % Lag range (-T/2)*dt, (T/2-1)*dt in ms
t([seg_size/2+1:seg_size,1:seg_size/2],2)=2*pi*real(ifft(fyx));  % Cumulant density, PBMB (5.9)
t([seg_size/2+1:seg_size,1:seg_size/2],3)=rhoyx;                 % rhoyx, JIN (2.8)

% Estimate variance of cumulant density and rho estimates
var_fac=4*pi*pi/(seg_size*samp_tot);                           % Factor (2pi/T)(2pi/R)
q_var=var_fac*2*sum(fxx(1:seg_size/2+1).*fyy(1:seg_size/2+1)); % PBMB (6.10)
rho_var=1/samp_tot;                                            % JIN (2.40)

cl.type=0;                   % Analysis type
cl.seg_size=seg_size;        % T
cl.seg_tot=seg_tot;          % L
cl.seg_tot_var=seg_tot;      % Effective no of segments, L'
cl.samp_tot=seg_tot*seg_size;% R=LT
cl.samp_rate=samp_rate;      % Sampling rate.
cl.dt=deltat;                % Delta t
cl.df=deltaf;                % Delta f
cl.f_c95=0.8512*sqrt(1/cl.seg_tot_var);  % 95% Confidence limit for spectral estimates, PBMB (6.2)
                                         % Confidence interval in plot is TWICE this value
cl.ch_c95=1-0.05^(1/(cl.seg_tot_var-1)); % 95% Confidence limit for coherence, PBMB (6.6)
cl.q_c95=1.96*sqrt(q_var);               % 95% Confidence limits for cumulant, PBMB (6.11)
cl.rho_c95=1.96*sqrt(rho_var);           % 95% Confidence limits for rho, JIN (2.41)
cl.col_R20=10;               % column in f matrix containing zero lag coherence component
cl.col_R2p=11;               % column in f matrix containing forward (positive lag) coherence component
cl.col_R2n=12;               % column in f matrix containing reverse (negative lag) coherence component
cl.col_rho=3;                % column in t matrix containing rho estimate
cl.R2=R2_rho;                % R2 value, R^2_yx
cl.R2_0=R2_rho_0;            % Component of R2 at zero lag, R^2_yx;0
cl.R2_p=R2_rho_p;            % Component of R2 in forward direction, positive lag component, R^2_yx;+
cl.R2_n=R2_rho_n;            % Component of R2 in reverse direction, negative lag component, R^2_yx;-
cl.R2_ch=R2_ch;              % R2 estimated from coherence between original processes
cl.R2_chw=R2_chw;            % R2 estimated from coherence between MMSE whitened processes
cl.R2_fprime_0=R2_fprime_0;  % Component of R2 at zero lag estimated from f'
cl.R2_fprime_p=R2_fprime_p;  % Component of R2 at positive lag estimated from f'
cl.R2_fprime_n=R2_fprime_n;  % Component of R2 at negative lag estimated from f'
cl.N2=0;        % N2, No of events in ch2, for point process data, zero here
cl.N1=0;        % N1, No of events in ch1, for point process data, zero here
cl.P2=0;        % P2, mean intensity of ch2, for point process data, zero here
cl.P1=0;        % P1, mean intensity of ch1, for point process data, zero here
cl.opt_str='';  % Options string, empty here
cl.what='';     % Field for plot label

% Display messages
% disp(['Segments: ',num2str(seg_tot),', Segment length: ',num2str(seg_size/samp_rate),' sec,  Resolution: ',num2str(cl.df),' Hz.']);
% disp(['  R2 (0, p, n): ',num2str(cl.R2),' (',num2str([cl.R2_0 cl.R2_p cl.R2_n]),').'])
