function [] = wrapper_Fig2_SymSNR(C,NCV,NC)
[bfilt,afilt] = buildfilter();
% NCvec = linspace(0,2,NC);
NCvec = [0 0.75 3];
NCvec = sqrt(NCvec); % variance
nc_col_sc = [1 0.9 0.8];
for i = 1:NC
    NCtits{i} = ['SigLeak = ' num2str(NCvec(i))];
end

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
    if ncov ==1
        bstrp = 0; %1;
    else
        bstrp = 0;
    end
    data = X;
    plotfig =1;
    linestyle = lsstyles{ncov};
    cmapn = cmap.*nc_col_sc(ncov);
    ncv = NCvec(ncov);
    
    %% Compute Signal Mixing
    randproc = randn(size(data.trial{1}));
    for i = 1:size(randproc,1)
        s = X.trial{1}(i,:);
        s = (s-mean(s))./std(s);
        n = ((NCvec(ncov)*1).*randproc(i,:));
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
    %     figure(1)
    %     plot(data.time{1},data.trial{1}+[0 2.5 5]');
    %     xlim([100 105])
    %     xlabel('Time'); ylabel('Amplitude'); grid on
    %     legend({'X1','X2','X3'})
    
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
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = computeNPD(data,1,8,1,bstrp);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm{1}; coh.ci = nscohspctrm{2};
    % NS Coh
    plotNSCoherence(coh,cmapn(1,:),Nsig,plotfig,linestyle,bstrp)
    %     plotNPD_zero(Hz,npdspctrm,data,cmap(1,:),plotfig,linestyle)
    % NPD
    plotNPD(Hz,npdspctrm,data,cmapn(3,:),plotfig,linestyle,bstrp)
    % NPDx1
    %     plotNPD(Hz,npdspctrmZ,data,cmapn(4,:),plotfig,linestyle)
    %     plotNPD(Hz,npdspctrmW,data,cmap(4,:),plotfig,linestyle)
    
    %% GRANGER
    figure(2)
    computeGranger(freq,cmapn(2,:),Nsig,plotfig,linestyle,0,bstrp)
    
    %% NPD CORR
    figure(3)
    plotNPDXCorr(lags,npdcrcv,data,[0 0 0],plotfig,linestyle)
    a =1;
end
SNRDB = 10*log10(snrbank(:,1));
SNRBPDB = 10*log10(snrbpbank(:,1));

a= 1;

% SNRDB:
%        Inf
%     1.2280
%    -4.7618

% SNRBPDB:
%        Inf
%     5.0973
%    -0.8948

% cfg = [];
% cfg.viewmode = 'butterfly';  % you can also specify 'butterfly'
% ft_databrowser(cfg, data);
