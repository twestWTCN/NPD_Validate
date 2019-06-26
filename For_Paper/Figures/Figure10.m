close all; clear
% addpath('C:\Users\Tim\Documents\Work\GIT\BrewerMap')
cmap = linspecer(4);

%% Simulation 10 - Combining ASNR and SigMix
rng(946022)
ncons = 4;
nreps = 25;
Nsig = 3;
[CMat,NCV] = makeRndGraphs(ncons,nreps,Nsig);

NC = 20;
% ASNR
SNRvec(1,:) = logspace(log10(0.01),log10(100),NC);
SNRvec(2,:) = repmat(log10(0.01),1,NC);
SNRvec(3,:) = repmat(log10(0.01),1,NC);
SNRvec = sqrt(SNRvec);
% SigMix
SigMixvec = linspace(0,1,NC);

% Nsamps = 150;

% Simulate data

for SNR = 1:size(SNRvec,2)
    for SMx = 1:size(SigMixvec,2)
    clear dataBank data
        parfor n = 1:nreps
            % Simulate Data
            cfg             = [];
            cfg.fsample     = 200;
            cfg.triallength = 500;
            cfg.ntrials     = 1;
            cfg.nsignal     = 3;
            cfg.method      = 'ar';
            cfg.params = CMat{2,n}; %selected two connections here!
            cfg.noisecov = NCV;
            data_raw = ft_connectivitysimulation(cfg);
            
            % Do the Signal Mixing
            data = data_raw;
            sigmix = repmat(SigMixvec(SMx)/(Nsig-1),Nsig,Nsig).*~eye(Nsig);
            sigmix = sigmix+eye(Nsig).*1; %(1-NCvec(ncov));
            data.trial{1} = (data.trial{1} -mean(data.trial{1},2))./std(data.trial{1},[],2);
            data.trial{1} = sigmix*data.trial{1};
            
            % NowDo the ASNR
            randproc = randn(size(data.trial{1}));
            for i = 1:size(randproc,1)
                s = data.trial{1}(i,:);
                s = (s-mean(s))./std(s);
                nr  =((SNRvec(i,SNR)*1).*randproc(i,:));
                y = s+nr;
                snr = var(s)/var(nr);
                disp(snr)
                %                 snrbank(SMx,n,i) = snr;
                %         y = (y-mean(y))./std(y);
                data.trial{1}(i,:) = y;
            end
            dataBank{n} = data;
            %              dataBank{SNR,SMx,n} = data;
            disp([SNR SMx n])
        end
    mkdir([cd '\benchmark'])
    save([cd '\benchmark\simdata_10_' sprintf('%.0f_%.0f',[SNR,SMx])],'dataBank')
    end
end

% Now test for recovery with dFC metrics
load([cd '\benchmark\10_NPG_CI_NPD_CI'],'NPG_ci','NPD_ci')

for SNR = 1:size(SNRvec,2)
    for SMx = 1:size(SigMixvec,2)
        bstrap = 0;
        load([cd '\benchmark\10_NPG_CI_NPD_CI'],'NPG_ci','NPD_ci')
        load([cd '\benchmark\simdata_10_' sprintf('%.0f_%.0f',[SNR,SMx])],'dataBank')
        parfor n = 1:nreps
            
            TrueCMat = CMat{2,n};
            Z = sum(TrueCMat,3);
            Z(Z==0.5) = 0;
            Z(Z==0.3) = 1;
            dataN = dataBank{n};
            %             dataN = dataBank{SNR,SMx,n};
            freq = computeSpectra(dataN,[0 0 0],3,0,'-',-1,(2^9)/dataN.fsample);
            [Hz granger grangerft] = computeGranger(freq,cmap(2,:),3,0,'--',1,bstrap);
            [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = computeNPD(dataN,1,9,1,bstrap);
%             plotNPD(Hz,npdspctrm,dataN,cmap(3,:),1,':',1)
            %             if bstrap == 1
            %                 NPG_ci{SNR,SMx}= granger{2,2};
            %                 NPD_ci{SNR,SMx} = npdspctrm{2,2};
            %                 bstrap = 0;
            %                 save([cd '\benchmark\10_NPG_CI_NPD_CI'],'NPG_ci','NPD_ci')
            %             else
            %                 load([cd '\benchmark\10_NPG_CI_NPD_CI'],'NPG_ci','NPD_ci')
            %             end
            
            % Now estimate
            NPG = granger{1,2};
            A = squeeze(sum((NPG>NPG_ci{4,4}),3));
            crit = ceil(size(NPG,3).*0.3);
            Ac = A>crit;
            NPGScore(SNR,SMx,n) = matrixScore(Ac,Z);
            
            NPD = npdspctrm{1,2};
            B = squeeze(sum((NPD>NPD_ci{4,4}),3));
            crit = ceil(size(NPD,3).*0.3);
            Bc = B>crit;
            NPDScore(SNR,SMx,n) = matrixScore(Bc,Z);
            %             NPDScore(dataLen,i,n) = sum((B(:)-Z(:)).^2);
            disp([SNR,SMx,n])
        end
    end
end

save([cd '\benchmark\10BenchMarks'],'NPGScore','NPDScore','SNRvec','SigMixvec')
load([cd '\benchmark\10BenchMarks'],'NPGScore','NPDScore','SNRvec','SigMixvec')

subplot(1,2,1);
[cma1,b1,cf] = contourf(log(SNRvec(1,:)),SigMixvec,mean(NPGScore,3),0:10:100);
clabel(cma1,b1)
title('NPG'); caxis([0 100])
xlabel('ASNR'); ylabel('Mixing (% Shared Variance)')
cmap = brewermap(10,'RdBu');
colormap(cmap)

subplot(1,2,2);
[cma2,b2,cf] = contourf(log(SNRvec(1,:)),SigMixvec,mean(NPDScore,3),0:10:100);
clabel(cma2,b2)
title('NPD'); caxis([0 100])
xlabel('ASNR'); ylabel('Mixing')
cmap = brewermap(10,'RdBu');
colormap(cmap)

set(gcf,'Position',[404         577        1109         419])

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
