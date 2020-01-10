function [] = wrapper_Fig5C_contSigmMix(C,NCV,NC)
% Sweep across signal mixings
NCvec = linspace(0,2,NC);
nc_col_sc = [1 0.9 0.8];

cmap = linspecer(4);
cmapn = cmap.*nc_col_sc(1);

% Run the MVAR Simulation (fieldtrip implementation uses BSMART)
rng(5231)
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
segOrd = 8; % 2^n length of segment used for FFT

for ncov = 1:NC
    disp(ncov)
    if ncov == 1
%         perm = 1;
%         permtype = 2;
        % Can switch off if you have run once to save the results see
        % L75-76
        perm = 0;
        permtype = 0; 
    else
        perm = 0;
        permtype = 0;
    end
    
        % Plotting parameters
    linestyle = '--';
    plotfig =0;

    %% Compute Signal Mixing
    data = X;
    sigmix = repmat(NCvec(ncov)/(Nsig-1),Nsig,Nsig).*~eye(Nsig);
    sigmix = sigmix+eye(Nsig).*1; %(1-NCvec(ncov));
    data.trial{1} = (data.trial{1} -mean(data.trial{1},2))./std(data.trial{1},[],2);
    data.trial{1} = sigmix*data.trial{1};
    dumdata = randn(size(data.trial{1}));
    shvar(:,:,ncov) = corr((sigmix*dumdata)');
    
    %% Compute Power with Fieldtrip (used for NPD and NPG computations)
    datalength = (2^segOrd)./cfg.fsample;
    freq = computeSpectra(data,[0 0 0],Nsig,plotfig,linestyle,1,datalength);
    pow = mean(abs(squeeze(freq.fourierspctrm(:,1,:))),1);
    npPow(ncov) = max(pow(freq.freq>42 & freq.freq<62));
    
    %% Compute NPD with Neurospec
    coh = [];
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = ft_computeNPD(freq,cfg.fsample,1,segOrd,perm,permtype);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm{1};
    nscoh(ncov) = max(nscohspctrm{1}(1,2,Hz>42 & Hz<62));
    npd(ncov) = max(npdspctrm{1,3}(1,2,Hz>42 & Hz<62));
    npdCi(ncov) = mean(npdspctrm{2,2}(1,2,:));

    %% Compute NPG
    [Hz granger grangerft] = computeGranger(freq,0,perm,permtype)
    npGC(ncov) = max(granger{1,3}(1,2,grangerft.freq>42 & grangerft.freq<62));
    npGCci(ncov) = mean(granger{2,3}(1,2,:));

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
plot(x,repmat(npdCi(1),1,size(x,2)),'LineStyle','--','color',cmapn(3,:))
plot(x,-repmat(npdCi(1),1,size(x,2)),'LineStyle','--','color',cmapn(3,:))

scatter(x,npGC,40,cmapn(2,:),'filled')
plot(x,repmat(npGCci(1),1,size(x,2)),'LineStyle','--','color',cmapn(2,:))
plot(x,-repmat(npGCci(1),1,size(x,2)),'LineStyle','--','color',cmapn(2,:))

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
