function [] = mvarconsim_npdver_F8(data,NC)

% NCvec = linspace(0,1,NC);
NCvec = [0 0.05 0.1];
nc_col_sc = [1 0.9 0.8];
Nsig = 2;
for i = 1:NC
    NCtits{i} = ['SigLeak = ' num2str(NCvec(i))];
end

cmap = linspecer(4);
lsstyles = {'-','-.',':'};
for ncov = 1:NC
    cmapn = cmap.*nc_col_sc(ncov);
    linestyle = lsstyles{ncov};
    
    %% Compute Signal Mixing
    sigmix = repmat(NCvec(ncov)/(Nsig-1),Nsig,Nsig).*~eye(Nsig);
    sigmix = sigmix+eye(Nsig).*(1-NCvec(ncov));
    data.trial{1} = sigmix*data.trial{1};
    data.trial{1} = (data.trial{1} -mean(data.trial{1},2))./std(data.trial{1},[],2);
    %     if ncov == NC
    %         plotfig =1;
    %     else
    plotfig =1;
    %     end
    %% Plot Example Trace
    figure(1)
    plot(data.time{1},data.trial{1}+[0 2.5]');
    xlim([100 105])
    xlabel('Time'); ylabel('Amplitude'); grid on
    legend({'X1','X2','X3'})
    
    %% Power
    figure(2)
    freq = computeSpectra(data,[0 0 0],Nsig,0,linestyle,3);
    
    %% Coherence
    %         figure(2)
    %    computeCoherence(freq,cmap(1,:),Nsig,plotfig,linestyle)
    
    %     %% WPLI
    %     computeWPLI(freq,cmap,Nsig,linestyle)
    
    %% NPD
    figure(2)
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = computeNPD(data,1,6);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm;
    % NS Coh
    computeNSCoherence(coh,cmapn(1,:),Nsig,0,linestyle)
    %     plotNPD_zero(Hz,npdspctrm,data,cmap(1,:),plotfig,linestyle)
    % NPD
    plotNPD(Hz,npdspctrm,data,cmapn(3,:),0,linestyle)
    % NPDx1
    %     plotNPD(Hz,npdspctrmZ,data,cmapn(4,:),plotfig,linestyle)
    %     plotNPD(Hz,npdspctrmW,data,cmap(4,:),plotfig,linestyle)
    
    %% GRANGER
    figure(2)
    [gHz granger] = computeGranger(freq,cmapn(2,:),Nsig,0,linestyle);
    
    plotSTN_M2_data(Hz,npdspctrm,data,cmapn(3,:),1,linestyle)
    plotSTN_M2_data(gHz,granger,data,cmapn(2,:),1,linestyle)
    
    
    
    %% NPD CORR
    %     figure(3)
    %     plotNPDXCorr(lags,npdcrcv,data,[0 0 0],plotfig,linestyle)

    
    a =1;
end


% cfg = [];
% cfg.viewmode = 'butterfly';  % you can also specify 'butterfly'
% ft_databrowser(cfg, data);
