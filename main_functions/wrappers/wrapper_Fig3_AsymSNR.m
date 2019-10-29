function [] = wrapper_Fig3_AsymSNR(C,NCV,NC)
% Asymetric SNRs
% NCvec = linspace(0,2,NC);
NCvec = [0.2  1      10
         0.2  0.2   0.2
         0.1  0.1   0.1];
NCvec = sqrt(NCvec);
nc_col_sc = [1 0.9 0.8];
for i = 1:NC
    NCtits{i} = ['SigLeak = ' num2str(NCvec(i))];
end
rng(8281)
cmap = linspecer(4);
lsstyles = {'-','-.',':'};
Nsig = size(C,1);
cfg             = [];
cfg.ntrials     = 1;
cfg.triallength = 250;
cfg.fsample     = 200;
cfg.nsignal     = Nsig;
cfg.method      = 'ar';
cfg.params = C;
cfg.noisecov = NCV;
X              = ft_connectivitysimulation(cfg);


for ncov = 1:NC
    if ncov == 1
        bstrp = 0;
    else
        bstrp = 0;
    end
    plotfig =1;
    linestyle = lsstyles{ncov};
    cmapn = cmap.*nc_col_sc(ncov);
    ncv = NCvec(ncov);
    data = X;
    %% Compute Signal Mixing
    randproc = randn(size(data.trial{1}));
    for i = 1:size(randproc,1)
        s = data.trial{1}(i,:);
        s = (s-mean(s))./std(s);
        n = ((NCvec(i,ncov)*1).*randproc(i,:));
        %         si = filter(bfilt,afilt,s);
        %         ni = filter(bfilt,afilt,n);
        y = s + n;
        snr = var(s)/var(n);
        snrbank(ncov,i) = snr;
        %         snrbp = var(si)/var(ni);
        snrbp = computeBandLimSNR(s,n,[45 55],data);
        snrbpbank(ncov,i) = snrbp;
        data.trial{1}(i,:) = y;

    end
    %% Plot Example Trace
    figure(1)
    plot(data.time{1},data.trial{1}+[0 2.5 5]');
    xlim([100 105])
    xlabel('Time'); ylabel('Amplitude'); grid on
    legend({'X1','X2','X3'})
    
    %% Power
    figure(2)
    freq = computeSpectra(data,[0 0 0],Nsig,plotfig,linestyle);
    
    %% Coherence
    %     figure(2)
    %     computeCoherence(freq,cmap(1,:),Nsig,plotfig,linestyle)
    
    %     %% WPLI
    %     computeWPLI(freq,cmap,Nsig,linestyle)
    
    %% NPD
    figure(2)
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = computeNPD(data,1);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm{1}; coh.ci = nscohspctrm{2};
    % NS Coh
    %     computeNSCoherence(coh,cmapn(1,:),Nsig,plotfig,linestyle)
    %     plotNPD_zero(Hz,npdspctrm,data,cmap(1,:),plotfig,linestyle)
    % NPD
%     plotNPD(Hz,npdspctrm,data,cmapn(3,:),plotfig,linestyle,bstrp)
    % NPDx1
    %     plotNPD(Hz,npdspctrmZ,data,cmapn(4,:),plotfig,linestyle)
    %     plotNPD(Hz,npdspctrmW,data,cmap(4,:),plotfig,linestyle)
    plotDiffNPD(Hz,npdspctrm,data,cmapn(3,:),plotfig,linestyle,1)
    %% GRANGER
    figure(2)
    computeGranger(freq,cmapn(2,:),Nsig,plotfig,linestyle,0,bstrp)
    
    %% NPD CORR
    %     figure(3)
    %     plotNPDXCorr(lags,npdcrcv,data,[0 0 0],plotfig,linestyle)
    a =1;
end

SNRDB = 10*log10(snrbpbank(:,1));
SNRBASE = 10*log10(snrbpbank(:,2));
SRAT = SNRDB-SNRBASE;

a= 1;


% SNRBASE =
% 
%    13.1839
%    13.4301
%    13.2043
% SNRDB =
% 
%    13.2253
%     6.1606
%    -3.8209
% SRAT =
% 
%     0.0414
%    -7.2695
%   -17.0252
% 

% cfg = [];
% cfg.viewmode = 'butterfly';  % you can also specify 'butterfly'
% ft_databrowser(cfg, data);
