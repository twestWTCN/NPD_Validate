% NPD_Validate_AddPaths()
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools')
clear; close all

%% This script will reproduce Figure 7 - Investigating the role of combined
% data confounds (instantaneous mixing and asymmetric signal-to-noise ratio;
% SNR_XY) upon the accuracy of connectivity estimation when using
% non-parametric directionality, and non-parametric multivariate Granger
% causality.

%% Simulation 10 - Combining ASNR and SigMix
% for permrun = 1:2
fresh = 1; % (1) Run permutation tests; (0) Load precomputed tables.
for permrun = 1:2 % (1) FFT Shuffle; (2) Phase Randomize
    
    % The length is fixed so use previously computed values
    %     load([cd '\benchmark\9A_NPG_CI_NPD_CI_' num2str(permtype)],'NPG_ci','NPD_ci')
    %     % This is 2^8 data length at 200s availability
    %     NPD_ci = NPD_ci{10}; NPG_ci = NPG_ci{10};
    %
    for Ncon = 1:2
        rng(946022)
        ncons = 3;
        %         Ncon = Ncon;
        nreps =24;
        Nsig = 3;
        [CMat,NCV] = makeRndGraphs(ncons,nreps,Nsig);
        
        NC = 13; % Needs to be Odd in order to get a scale
        % ASNR
        SNRvec = repmat(0.001,Nsig,NC);
        SNRvec(1,:) = logspace(log10(1000),log10(0.01),NC);
%         SNRvec(1,1:floor(NC/2)) = logspace(log10(1000),log10(0.01),floor(NC/2));
%         SNRvec(2,ceil(NC/2)+1:end) = logspace(log10(0.01),log10(1000),floor(NC/2));
        SNRvec = sqrt(SNRvec);
        % SigMix
        SigMixvec = linspace(0,1,13);
        % Nsamps = 150;
        
%         Simulate data
        if fresh == 1
            for SNR = 1:size(SNRvec,2)
                for SMx = 1:size(SigMixvec,2)
                    clear dataBank data
                    parfor n = 1:nreps
                        % Simulate Data
                        cfg             = [];
                        cfg.fsample     = 200;
                        cfg.triallength = 200;
                        cfg.ntrials     = 1;
                        cfg.nsignal     = 3;
                        cfg.method      = 'ar';
                        cfg.params = CMat{Ncon,n}; %selected Ncon connections here!
                        cfg.noisecov = NCV;
                        data_raw = ft_connectivitysimulation(cfg);
                        
                        % Find Connectivity Matrix (to bias one source)
                        TrueCMat = CMat{Ncon,n};
                        Z = sum(TrueCMat,3);
                        Z(Z==0.5) = 0;
                        Z(Z~=0) = 1;
                        Z = Z>0;
                        [src,tar] = find(Z);
                        AOrd = [src(1) tar(1) setdiff(1:3,[src(1) tar(1)])]; % This sets the order of SNR changes
                        
                        % Do the Signal Mixing
                        data = data_raw;
                        sigmix = repmat(SigMixvec(SMx)/(Nsig-1),Nsig,Nsig).*~eye(Nsig);
                        sigmix = sigmix+eye(Nsig).*1; %(1-NCvec(ncov));
                        data.trial{1} = (data.trial{1} -mean(data.trial{1},2))./std(data.trial{1},[],2);
                        data.trial{1} = sigmix*data.trial{1};
                        
                        % NowDo the ASNR
                        randproc = randn(size(data.trial{1}));
                        snrbp = [];
                        for i = AOrd
                            s = data.trial{1}(i,:);
                            s = (s-mean(s))./std(s);
                            nr  =((SNRvec(i,SNR)*1).*randproc(i,:));
                            y = s+nr;
                            snr = var(s)/var(nr);
                            disp(snr)
                            snrbp(i) = computeBandLimSNR(s,nr,[45 55],data);
                            data.trial{1}(i,:) = y;
                        end
                        SNRDB = 10*log10(snrbp);
                        
                        deltaDB = SNRDB(1) - min(SNRDB(2:3));
                        deltaDBBank(SNR,SMx,n) = deltaDB;
                        %                         deltaDBBank(n) = deltaDB;
                        dataBank{n} = data;
                        %              dataBank{SNR,SMx,n} = data;
                        disp([SNR SMx n])
                    end
                    mkdir([cd '\benchmark'])
                    save([cd '\benchmark\simdata_10_' sprintf('%.0f_%.0f_%.0f',[SNR,SMx,Ncon])],'dataBank','deltaDBBank','CMat')
                end
            end
        end
        
        segOrd = 8; % 2^n length of segment used for FFT
        
        % Compute CI Threshold
        perm = 0;
        if perm == 1
            permtype = permrun;
            plotfig = 0;
            linestyle = '-';
            load([cd '\benchmark\simdata_10_' sprintf('%.0f_%.0f_%.0f',[1,1,Ncon])],'dataBank','deltaDBBank','CMat')
            dataN = dataBank{1};
            datalength = (2^segOrd)./dataN.fsample;
            freq = computeSpectra(dataN,[0 0 0],Nsig,plotfig,linestyle,-1,datalength);
            [Hz granger grangerft] = computeGranger(freq,1,perm,permtype)
            [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = ft_computeNPD(freq,dataN.fsample,1,segOrd,perm,permtype);
            
            NPG_ci{1}= granger{2,2};
            NPD_ci{1} = npdspctrm{2,2};
            perm = 0;
            save([cd '\benchmark\10_NPG_CI_NPD_CI_' num2str(permrun) '_NCon' num2str(Ncon)],'NPG_ci','NPD_ci')
        else
            permtype = permrun;
            load([cd '\benchmark\10_NPG_CI_NPD_CI_' num2str(permrun) '_NCon' num2str(Ncon)],'NPG_ci','NPD_ci')
        end
        
        
        
        % Now test for recovery with dFC metrics
        % load([cd '\benchmark\10_NPG_CI_NPD_CI_' num2str(Ncon)],'NPG_ci','NPD_ci')
        NPDScore = []; NPGScore = []; SNRRec = [];
        for SNR = 1:size(SNRvec,2)
            for SMx = 1:size(SigMixvec,2)
                load([cd '\benchmark\simdata_10_' sprintf('%.0f_%.0f_%.0f',[SNR,SMx,Ncon])],'dataBank','deltaDBBank','CMat')
                load([cd '\benchmark\10_NPG_CI_NPD_CI_' num2str(permrun) '_NCon' num2str(Ncon)],'NPG_ci','NPD_ci')
                parfor n = 1:nreps % if the CI is precomputed
                    %                 for n = 1:nreps
                    TrueCMat = CMat{Ncon,n};
                    Z = sum(TrueCMat,3);
                    Z(Z==0.5) = 0;
                    Z(Z~=0) = 1;
                    Z = Z>0;
                    dataN = dataBank{n};
                    plotfig = 0;
                    linestyle = '-';
                    %             dataN = dataBank{SNR,SMx,n};
                    
                    % Compute Spectra
                    datalength = (2^segOrd)./dataN.fsample;
                    freq = computeSpectra(dataN,[0 0 0],Nsig,plotfig,linestyle,-1,datalength);
                    % Compute connectivity
                    [Hz granger grangerft] = computeGranger(freq,1,0,permtype)
                    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = ft_computeNPD(freq,dataN.fsample,1,segOrd,0,permtype);
                    %             plotNPD(Hz,npdspctrm,dataN,'r',1,':',1)
                    
                    % Now estimateNPG_ci
                    NPG = granger{1,2};
                    A = squeeze(sum((NPG>NPG_ci{1}),3));
                    crit = ceil(size(NPG,3).*0.1);
                    Ac = A>crit;
                    NPGScore(SNR,SMx,n) = matrixScore(Ac,Z);
                    
                    NPD = npdspctrm{1,2};
                    B = squeeze(sum((NPD>NPD_ci{1}),3));
                    crit = ceil(size(NPD,3).*0.1);
                    Bc = B>crit;
                    NPDScore(SNR,SMx,n) = matrixScore(Bc,Z);
                    deltaDBScore(SNR,SMx,n) = deltaDBBank(SNR,SMx,n)
                    %             NPDScore(dataLen,i,n) = sum((B(:)-Z(:)).^2);
                    disp([SNR,SMx,n])
                end
            end
        end
        
        save([cd '\benchmark\10BenchMarks_' num2str(Ncon) '_' num2str(permtype)],'NPGScore','NPDScore','deltaDBScore','SNRvec','SigMixvec')
        load([cd '\benchmark\10BenchMarks_' num2str(Ncon) '_' num2str(permtype)],'NPGScore','NPDScore','deltaDBScore','SNRvec','SigMixvec')
        
    end
    
end

for X = 1:2
    for Nc = 1:2
        
        figure(X*10 + Nc)
        
        load([cd '\benchmark\10BenchMarks_' num2str(Nc) '_' num2str(X)],'NPGScore','NPDScore','deltaDBScore','SNRvec','SigMixvec')
        
        meandeltaDB = nanmean(squeeze(deltaDBScore(:,1,:)),2);
%         meandeltaDB = mean(mean(deltaDBBank,3),2);
        subplot(1,2,1);
        [cma1,b1,cf] = contourf(meandeltaDB,SigMixvec,mean(NPGScore,3)',-10:10:100);
        clabel(cma1,b1)
        title('NPG');
        % caxis([-10 75])
        caxis([-10 100])
        
        xlabel('ASNR'); ylabel('Mixing (% Shared Variance)')
        cmap = brewermap(21,'RdBu');
        colormap(cmap)
        xlim([-50 0])
        
        subplot(1,2,2);
        [cma2,b2,cf] = contourf(meandeltaDB,SigMixvec,mean(NPDScore,3)',-10:10:100);
        clabel(cma2,b2)
        title('NPD');
        caxis([-10 100])
        % caxis([-10 75])
        xlabel('ASNR'); ylabel('Mixing')
        colormap(cmap)
        colorbar
        set(gcf,'Position',[404         577        1109         419])
%         xlim([-50 0])
    end
end
% Bc = [0 0 1;
%       0 0 0;
%       0 0 0];
% Z  = [0 0 1;
%       1 0 0;
%       1 1 1];
%
%  matrixScore(Bc,Z)
%
%
