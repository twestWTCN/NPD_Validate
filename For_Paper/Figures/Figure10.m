close all; clear
% addpath('C:\Users\Tim\Documents\Work\GIT\BrewerMap')
cmap = linspecer(4);

%% Simulation 10 - Combining ASNR and SigMix
rng(946022)
ncons = 4;
Ncon = 1;
nreps = 24;
Nsig = 3;
[CMat,NCV] = makeRndGraphs(ncons,nreps,Nsig);

NC = 12;
% ASNR
SNRvec(1,:) = logspace(log10(0.001),log10(1000),NC);
SNRvec(2,:) = repmat(0.001,1,NC);
SNRvec(3,:) = repmat(0.001,1,NC);
% SNRvec(4,:) = repmat(0.001,1,NC);
% SNRvec = sqrt(SNRvec);
% SigMix
SigMixvec = linspace(0,1,NC);

% Nsamps = 150;

% % Simulate data

for SNR = 1:size(SNRvec,2)
    for SMx = 1:size(SigMixvec,2)
    clear dataBank data
        parfor n = 1:nreps
            % Simulate Data
            cfg             = [];
            cfg.fsample     = 200;
            cfg.triallength = 100;
            cfg.ntrials     = 1;
            cfg.nsignal     = 3;
            cfg.method      = 'ar';
            cfg.params = CMat{Ncon,n}; %selected Ncon connections here!
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
    save([cd '\benchmark\simdata_10_' sprintf('%.0f_%.0f_%.0f',[SNR,SMx,Ncon])],'dataBank')
    end
end

% Now test for recovery with dFC metrics
% load([cd '\benchmark\10_NPG_CI_NPD_CI_' num2str(Ncon)],'NPG_ci','NPD_ci')

for SNR = 1:size(SNRvec,2)
    for SMx = 1:size(SigMixvec,2)
        bstrap = 0;
%         load([cd '\benchmark\10_NPG_CI_NPD_CI_' num2str(Ncon)],'NPG_ci','NPD_ci')
        load([cd '\benchmark\simdata_10_' sprintf('%.0f_%.0f_%.0f',[SNR,SMx,Ncon])],'dataBank')
        parfor n = 1:nreps
            
            TrueCMat = CMat{Ncon,n};
            Z = sum(TrueCMat,3);
            Z(Z==0.5) = 0;
            Z(Z~=0) = 1;
            dataN = dataBank{n};
            %             dataN = dataBank{SNR,SMx,n};
            freq = computeSpectra(dataN,[0 0 0],3,0,'-',-1,(2^8)/dataN.fsample);
            [Hz granger grangerft] = computeGranger(freq,cmap(2,:),3,0,'--',1,bstrap);
            [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = computeNPD(dataN,1,8,1,bstrap);
%             plotNPD(Hz,npdspctrm,dataN,cmap(3,:),1,':',1)
            %             if bstrap == 1
            %                 NPG_ci{SNR,SMx}= granger{2,2};
            %                 NPD_ci{SNR,SMx} = npdspctrm{2,2};
            %                 bstrap = 0;
            %                 save([cd '\benchmark\10_NPG_CI_NPD_CI_' num2str(Ncon)],'NPG_ci','NPD_ci')
            %             else
            %                 load([cd '\benchmark\10_NPG_CI_NPD_CI_' num2str(Ncon)],'NPG_ci','NPD_ci')
            %             end
            %
            % Now estimate
            NPG = granger{1,2};
            A = squeeze(sum((abs(NPG)>0.05),3));
            crit = ceil(size(NPG,3).*0.15);
            Ac = A>crit;
            NPGScore(SNR,SMx,n) = matrixScore(Ac,Z);
            
            NPD = npdspctrm{1,2};
            B = squeeze(sum((NPD>0.05),3));
            crit = ceil(size(NPD,3).*0.15);
            Bc = B>crit;
            
            NPD_z = npdspctrmZ{1,2};
            B = squeeze(sum((NPD_z>0.05),3));
            crit = ceil(size(NPD_z,3).*0.15);
            Bc_z = B>crit;
            
            NPD_w = npdspctrmW{1,2};
            B = squeeze(sum((NPD_w>0.05),3));
            crit = ceil(size(NPD_w,3).*0.15);
            Bc_w = B>crit;
 
%             NPD_q = npdspctrmQ{1,3};
%             B = squeeze(sum((NPD_q>0.05),3));
%             crit = ceil(size(NPD_q,3).*0.15);
%             Bc_q = B>crit;
            
            NPDScore(SNR,SMx,n) = matrixScore(Bc,Z);
            
            %             NPDScore(dataLen,i,n) = sum((B(:)-Z(:)).^2);
            disp([SNR,SMx,n])
        end
    end
end

save([cd '\benchmark\10BenchMarks_' num2str(Ncon)],'NPGScore','NPDScore','SNRvec','SigMixvec')
load([cd '\benchmark\10BenchMarks_' num2str(Ncon)],'NPGScore','NPDScore','SNRvec','SigMixvec')

subplot(1,2,1);
[cma1,b1,cf] = contourf(log(SNRvec(1,:)),SigMixvec,mean(NPGScore,3),-25:10:35);
clabel(cma1,b1)
title('NPG');
% caxis([-10 75])
caxis([-25 30])

xlabel('ASNR'); ylabel('Mixing (% Shared Variance)')
cmap = brewermap(21,'RdBu');
colormap(cmap)

subplot(1,2,2);
[cma2,b2,cf] = contourf(log(SNRvec(1,:)),SigMixvec,mean(NPDScore,3),-25:10:35);
clabel(cma2,b2)
title('NPD');
caxis([-25 30])
% caxis([-10 75])
xlabel('ASNR'); ylabel('Mixing')
colormap(cmap)
colorbar
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
