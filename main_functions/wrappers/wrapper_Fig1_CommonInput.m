function [] = wrapper_Fig1_CommonInput(C,NCV,NC,frstord,sndord,fname)

% NCvec = linspace(0,1,NC);
NCvec = [0 0.5 0.8];
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
data              = ft_connectivitysimulation(cfg);

for ncov = 1:NC
    if ncov == 1
        bstrp = 1;
    else
        bstrp = 0;
    end
    cmapn = cmap.*nc_col_sc(ncov);
    ncv = NCvec(ncov);
    
    linestyle = lsstyles{ncov};
    % Weighted mixture of signals
    %     sigmix = [
    %         0.7 0.1 0.1 0.1;
    %         0.1 0.7 0.1 0.1;
    %         0.1 0.1 0.7 0.1;
    %         0.1 0.1 0.1 0.7;
    %         ];
    %% Compute Signal Mixing
    sigmix = repmat(NCvec(ncov)/(Nsig-1),Nsig,Nsig).*~eye(Nsig);
    sigmix = sigmix+eye(Nsig).*(1-NCvec(ncov));
    
    data.trial{1} = sigmix*data.trial{1};
    
    
    %     if ncov == NC
    %         plotfig =1;
    %     else
    plotfig =1;
    %     end
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
    %         figure(2)
    %    computeCoherence(freq,cmap(1,:),Nsig,plotfig,linestyle)
    
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
    plotNPD(Hz,npdspctrmZ,data,cmapn(4,:),plotfig,linestyle,bstrp)
    %     plotNPD(Hz,npdspctrmW,data,cmap(4,:),plotfig,linestyle)
    
    %% GRANGER
    figure(2)
    computeGranger(freq,cmapn(2,:),Nsig,plotfig,linestyle,1,bstrp)
    
    %% NPD CORR
    %     figure(3)
    %     plotNPDXCorr(lags,npdcrcv,data,[0 0 0],plotfig,linestyle)
    
    a =1;
end
% save('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\precomp_CI_table\F1_CItab','npdCi','npGCci')
load('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\precomp_CI_table\F1A_CItab','npdCi','npGCci')


% cfg = [];
% cfg.viewmode = 'butterfly';  % you can also specify 'butterfly'
% ft_databrowser(cfg, data);
