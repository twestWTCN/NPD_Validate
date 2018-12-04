function [] = mvarconsim_npdver_F6(C,NCV,NC)

% NCvec = linspace(0,2,NC);
NCvec = [0 0.5 1.5];
nc_col_sc = [1 0.9 0.8];
% for i = 1:NC
%     NCtits{i} = ['SigLeak = ' num2str(NCvec(i))];
% end

cmap = linspecer(4);
lsstyles = {'-','-.',':'};
Nsig = size(C,1);
cfg             = [];
cfg.ntrials     = 1;
cfg.triallength = 500;
cfg.fsample     = 100;
cfg.nsignal     = Nsig;
cfg.method      = 'ar';
cfg.params = C;
cfg.noisecov = NCV;
data              = ft_connectivitysimulation(cfg);


for ncov = 1:NC
        plotfig =0;
    linestyle = lsstyles{ncov};
    cmapn = cmap.*nc_col_sc(ncov);
    ncv = NCvec(ncov);
    
    randproc = randn(size(data.trial{1}));
    for i = 2
        y = data.trial{1}(i,:);
        y = y+((NCvec(ncov)*std(y)).*randproc(i,:));
        y = (y-mean(y))./std(y);
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
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv npdcrcvZ npdcrcvW] = computeNPD(data,1);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm;
    % NS Coh
%     computeNSCoherence(coh,cmapn(1,:),Nsig,plotfig,linestyle)
    %     plotNPD_zero(Hz,npdspctrm,data,cmap(1,:),plotfig,linestyle)
    % NPD
    plotNPD(Hz,npdspctrm,data,cmapn(3,:),plotfig,linestyle)
    % NPDx1
    plotNPD(Hz,npdspctrmW,data,cmapn(4,:),plotfig,linestyle)
    %     plotNPD(Hz,npdspctrmW,data,cmap(4,:),plotfig,linestyle)
    
    %% GRANGER
    figure(2)
    computeGranger(freq,cmapn(2,:),Nsig,plotfig,linestyle)
    
    %% NPD CORR
        figure(3)
        plotNPDXCorr(lags,npdcrcv,data,cmapn(3,:),0,linestyle)
        plotNPDXCorr(lags,npdcrcvW,data,cmapn(4,:),0,linestyle)
    a =1;
end


% cfg = [];
% cfg.viewmode = 'butterfly';  % you can also specify 'butterfly'
% ft_databrowser(cfg, data);
