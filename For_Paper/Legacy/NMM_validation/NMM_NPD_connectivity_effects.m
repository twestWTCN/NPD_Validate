function [] = NMM_NPD_connectivity_effects(dataS,NC,frstord,sndord,fname)

% % NCvec = linspace(0,1,NC);
% NCvec = [0 0.5 0.8];
% nc_col_sc = [1 0.9 0.8];
% for i = 1:NC
%     NCtits{i} = ['SigLeak = ' num2str(NCvec(i))];
% end

cmapn = linspecer(4);
lsstyles = {'-','--','-.',':'};
for ncov = 1:NC
    linestyle = lsstyles{ncov};
    plotfig =1;
    data = dataS(ncov);
    Nsig = size(data.trial{1},1);
    %% Plot Example Trace
    figure(1*(10*ncov))
    plot(data.time{1},data.trial{1}+[0 2.5 5]');
    xlim([15 20])
    xlabel('Time'); ylabel('Amplitude'); grid on
    legend({'X1','X2','X3'})
    
    %% Power
    figure(2+(10*ncov))
    computeSpectra(data,[0 0 0],Nsig,plotfig,linestyle);
    freq = computeSpectra(data,[0 0 0],Nsig,0,linestyle,10);
    %% Coherence
    %         figure(2)
    %    computeCoherence(freq,cmap(1,:),Nsig,plotfig,linestyle)
    
    %     %% WPLI
    %     computeWPLI(freq,cmap,Nsig,linestyle)
    
    %% NPD
    figure(2+(10*ncov))
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = computeNPD(data,frstord,9);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm;
    % NS Coh
%     computeNSCoherence(coh,cmapn(1,:),Nsig,plotfig,linestyle)
        plotNPD_zero(Hz,npdspctrm,data,cmapn(1,:),1,linestyle)
    % NPD
    plotNPD(Hz,npdspctrm,data,cmapn(3,:),plotfig,linestyle)
    % NPDx1
    plotNPD(Hz,npdspctrmZ,data,cmapn(4,:),plotfig,linestyle)
    %     plotNPD(Hz,npdspctrmW,data,cmap(4,:),plotfig,linestyle)
    
    %% GRANGER
    figure(2+(10*ncov))
%     computeGranger(freq,cmapn(2,:),Nsig,plotfig,linestyle)
    
    %% NPD CORR
    figure(3+(10*ncov))
    plotNPDXCorr(lags,npdcrcv,data,[0 0 0],0,linestyle)
    
    
    a =1;
end


% cfg = [];
% cfg.viewmode = 'butterfly';  % you can also specify 'butterfly'
% ft_databrowser(cfg, data);
