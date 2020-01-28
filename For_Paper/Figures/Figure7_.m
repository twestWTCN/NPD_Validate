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
for permtype = 1:2 % (1) FFT Shuffle; (2) Phase Randomize
    
    % The length is fixed so use previously computed values
    load([cd '\benchmark\9A_NPG_CI_NPD_CI_' num2str(permtype)],'NPG_ci','NPD_ci')
    % This is 2^8 data length at 200s availability
    NPD_ci = NPD_ci{5}; NPG_ci = NPG_ci{5};
    
    for Ncon = 1:2
        rng(946022)
        ncons = 4;
        %         Ncon = Ncon;
        nreps = 24;
        Nsig = 3;
        [CMat,NCV] = makeRndGraphs(ncons,nreps,Nsig);
        
        NC = 12;
        % ASNR
        SNRvec(1,:) = logspace(log10(0.001),log10(1000),NC);
        SNRvec(2,:) = repmat(0.2,1,NC);
        SNRvec(3,:) = repmat(0.2,1,NC);
        % SNRvec(4,:) = repmat(0.001,1,NC);
        % SNRvec = sqrt(SNRvec);
        % SigMix
        SigMixvec = linspace(0,1,NC);
        % Nsamps = 150;
        
        %     Simulate data
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
                        AOrd = [src(1) setdiff(1:3,src(1))];
                        
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
                        deltaDB = SNRDB(tar(1))-mean(SNRDB(setdiff(1:Nsig,tar)));
                        deltaDBBank(n) = deltaDB;
                        dataBank{n} = data;
                        %              dataBank{SNR,SMx,n} = data;
                        disp([SNR SMx n])
                    end
                    mkdir([cd '\benchmark'])
                    save([cd '\benchmark\simdata_10_' sprintf('%.0f_%.0f_%.0f',[SNR,SMx,Ncon])],'dataBank','deltaDBBank','CMat')
                end
            end
        end
        
        
        % Now test for recovery with dFC metrics
        % load([cd '\benchmark\10_NPG_CI_NPD_CI_' num2str(Ncon)],'NPG_ci','NPD_ci')
        segOrd = 8; % 2^n length of segment used for FFT
        NPDScore = []; NPGScore = []; SNRRec = [];
        for SNR = 1:size(SNRvec,2)
            for SMx = 1:size(SigMixvec,2)
                load([cd '\benchmark\simdata_10_' sprintf('%.0f_%.0f_%.0f',[SNR,SMx,Ncon])],'dataBank','deltaDBBank','CMat')
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
                    A = squeeze(sum((NPG>NPG_ci),3));
                    crit = ceil(size(NPG,3).*0.1);
                    Ac = A>crit;
                    NPGScore(SNR,SMx,n) = matrixScore(Ac,Z);
                    
                    NPD = npdspctrm{1,2};
                    B = squeeze(sum((NPD>NPD_ci),3));
                    crit = ceil(size(NPD,3).*0.1);
                    Bc = B>crit;
                    NPDScore(SNR,SMx,n) = matrixScore(Bc,Z);
                    deltaDBScore(SNR,SMx,n) = deltaDBBank(n)
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
        
        meandeltaDB = nanmean(nanmean(deltaDBScore,3),2);
        subplot(1,2,1);
        [cma1,b1,cf] = contourf(meandeltaDB,SigMixvec,mean(NPGScore,3),-60:10:100);
        clabel(cma1,b1)
        title('NPG');
        % caxis([-10 75])
        caxis([-50 100])
        
        xlabel('ASNR'); ylabel('Mixing (% Shared Variance)')
        cmap = brewermap(21,'RdBu');
        colormap(cmap)
        
        subplot(1,2,2);
        [cma2,b2,cf] = contourf(meandeltaDB,SigMixvec,mean(NPDScore,3),-60:10:100);
        clabel(cma2,b2)
        title('NPD');
        caxis([-50 100])
        % caxis([-10 75])
        xlabel('ASNR'); ylabel('Mixing')
        colormap(cmap)
        colorbar
        set(gcf,'Position',[404         577        1109         419])
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
