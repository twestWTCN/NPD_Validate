function [] = wrapper_Fig5C_contSigmMix(C,NCV,NC)

NCvec = linspace(0,2,NC);
nc_col_sc = [1 0.9 0.8];
for i = 1:NC
    NCtits{i} = ['SigLeak = ' num2str(NCvec(i))];
end
rng(2321)

cmap = linspecer(4);
lsstyles = {'-','-.',':'};
cmapn = cmap.*nc_col_sc(1);
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
X = ft_connectivitysimulation(cfg);
for ncov = 1:NC
    if ncov == floor(NC/2)
        bstrp = 1;
    else
        bstrp = 0;
    end
    linestyle = '--';
    % Weighted mixture of signals
    %     sigmix = [
    %         0.7 0.1 0.1 0.1;
    %         0.1 0.7 0.1 0.1;
    %         0.1 0.1 0.7 0.1;
    %         0.1 0.1 0.1 0.7;
    %         ];
    %% Compute Signal Mixing
    data = X;
    sigmix = repmat(NCvec(ncov)/(Nsig-1),Nsig,Nsig).*~eye(Nsig);
    sigmix = sigmix+eye(Nsig).*1; %(1-NCvec(ncov));
    data.trial{1} = (data.trial{1} -mean(data.trial{1},2))./std(data.trial{1},[],2);
    data.trial{1} = sigmix*data.trial{1};
%     shvar(:,:,ncov) = corr(data.trial{1}');
    
    dumdata = randn(size(data.trial{1}));
    shvar(:,:,ncov) = corr((sigmix*dumdata)');
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
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = computeNPD(data,1,8,1,bstrp);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm{1};
    nscoh(ncov) = max(nscohspctrm{1}(1,2,Hz>42 & Hz<62));
    npd(ncov) = max(npdspctrm{1,3}(1,2,Hz>42 & Hz<62));
    npdCi(ncov) = mean(npdspctrm{2,2}(1,2,:));
    % NS Coh
    %     plotNSCoherence(coh,cmapn(1,:),Nsig,plotfig,linestyle)
    plotNPD_zero(Hz,npdspctrm,data,cmap(1,:),plotfig,linestyle)
    % NPD
    plotNPD(Hz,npdspctrm,data,cmapn(3,:),plotfig,linestyle)
    % NPDx1
    %     plotNPD(Hz,npdspctrmZ,data,cmapn(4,:),plotfig,linestyle)
    %     plotNPD(Hz,npdspctrmW,data,cmap(4,:),plotfig,linestyle)
    
    %% GRANGER
    [~, granger grangerft] = computeGranger(freq,cmapn(2,:),Nsig,plotfig,linestyle,1,bstrp);
    npGC(ncov) = max(granger{1,3}(1,2,grangerft.freq>42 & grangerft.freq<62));
    npGCci(ncov) = mean(granger{2,2}(1,2,:));
end
% save('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\precomp_CI_table\F5C_CItab','npdCi','npGCci')
load('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\precomp_CI_table\F5C_CItab','npdCi','npGCci')

for i = 1:size(shvar,3)
    p = shvar(:,:,i).*abs(eye(3)-1);
    x(i) =  sum(p(:))./sum(abs(p(:))>0);
end
a =1;
figure(3)
% scatter(NCvec,npPow,40,cmapn(1,:),'filled')
scatter(x,nscoh,50,cmapn(1,:),'Marker','+','LineWidth',3);
hold on
scatter(x,npd,40,cmapn(3,:),'filled')
plot(x,repmat(npdCi(floor(NC/2)),1,size(x,2)),'LineStyle','--','color',cmapn(3,:))
plot(x,-repmat(npdCi(floor(NC/2)),1,size(x,2)),'LineStyle','--','color',cmapn(3,:))

scatter(x,npGC,40,cmapn(2,:),'filled')
plot(x,repmat(npGCci(floor(NC/2)),1,size(x,2)),'LineStyle','--','color',cmapn(2,:))
plot(x,-repmat(npGCci(floor(NC/2)),1,size(x,2)),'LineStyle','--','color',cmapn(2,:))

grid on
xlabel('Shared Variance (correlation)');ylabel('FC Magnitude')
legend({'Coherence','NPD','Granger'})
title('FC vs Mixing')
% cfg = [];
% cfg.viewmode = 'butterfly';  % you can also specify 'butterfly'
% ft_databrowser(cfg, data);
figure(4)
% scatter(NCvec,npPow,40,cmapn(1,:),'filled')
scatter(NCvec,nscoh,50,cmapn(1,:),'Marker','+','LineWidth',3);
hold on
scatter(NCvec,npd,40,cmapn(3,:),'filled')
scatter(NCvec,npGC,40,cmapn(2,:),'filled')
grid on
xlabel('\lambda (mixing)');ylabel('FC Magnitude')
legend({'Coherence','NPD','Granger'})
title('FC vs Mixing')
