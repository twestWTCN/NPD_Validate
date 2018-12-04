function [] = mvarconsim_npdver_F2b(C,NCV,NC)

NCvec = linspace(0,0.5,NC);
nc_col_sc = [1 0.9 0.8];
for i = 1:NC
    NCtits{i} = ['SigLeak = ' num2str(NCvec(i))];
end

cmap = linspecer(4);
lsstyles = {'-','-.',':'};
cmapn = cmap.*nc_col_sc(1)
% ncv = NCvec(ncov);
Nsig = size(C,1);
cfg             = [];
cfg.ntrials     = 1;
cfg.triallength = 500;
cfg.fsample     = 200;
cfg.nsignal     = Nsig;
cfg.method      = 'ar';
cfg.params = C;
cfg.noisecov = NCV;
data              = ft_connectivitysimulation(cfg);
for ncov = 1:NC
    
    linestyle = '--';
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
    plotfig =0;
    %     end
    
    %% Power
    freq = computeSpectra(data,[0 0 0],Nsig,plotfig,linestyle);
        pow = mean(abs(squeeze(freq.fourierspctrm(:,1,:))),1);
    npPow(ncov) = max(pow(freq.freq>42 & freq.freq<62));
   
    %% NPD
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = computeNPD(data,1);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm;
    nscoh(ncov) = max(nscohspctrm(2,1,Hz>42 & Hz<62));
    npd(ncov) = max(npdspctrm{2}(2,1,Hz>42 & Hz<62));
    % NS Coh
    %     computeNSCoherence(coh,cmapn(1,:),Nsig,plotfig,linestyle)
    plotNPD_zero(Hz,npdspctrm,data,cmap(1,:),plotfig,linestyle)
    % NPD
    plotNPD(Hz,npdspctrm,data,cmapn(3,:),plotfig,linestyle)
    % NPDx1
    %     plotNPD(Hz,npdspctrmZ,data,cmapn(4,:),plotfig,linestyle)
    %     plotNPD(Hz,npdspctrmW,data,cmap(4,:),plotfig,linestyle)
    
    %% GRANGER
    granger = computeGranger(freq,cmapn(2,:),Nsig,plotfig,linestyle);
    npGC(ncov) = max(granger.grangerspctrm(1,2,granger.freq>42 & granger.freq<62));
end

a =1;
figure(1)
% scatter(NCvec,npPow,40,cmapn(1,:),'filled')
scatter(NCvec,nscoh,40,cmapn(1,:),'filled')
hold on
scatter(NCvec,npd,40,cmapn(3,:),'filled')
scatter(NCvec,npGC,40,cmapn(2,:),'filled')
grid on
xlabel('\lambda (mixing)');ylabel('FC Magnitude')
legend({'Coherence','NPD','Granger'})
title('FC vs Mixing')
% cfg = [];
% cfg.viewmode = 'butterfly';  % you can also specify 'butterfly'
% ft_databrowser(cfg, data);
