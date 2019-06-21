close all; clear
addpath('C:\Users\Tim\Documents\Work\GIT\BrewerMap')
cmap = linspecer(4);

%% Simulation 10 - Combining ASNR and SigMix
rng(946022)
ncons = 4;
nreps = 10;
Nsig = 3;
[CMat,NCV] = makeRndGraphs(ncons,nreps,Nsig);

NC = 5;
% ASNR 
SNRvec(1,:) = logspace(log10(0.01),log10(1000),NC);
SNRvec(2,:) = ones(1,NC);
SNRvec(3,:) = ones(1,NC);

% SigMix
SigMixvec = linspace(0,2,NC);

% Nsamps = 150;

% Simulate data

% for SNR = 1:size(SNRvec,2)
%     clear data
%     for SMx = 1:size(SigMixvec,2)
%         parfor n = 1:nreps
%             % Simulate Data
%             cfg             = [];
%             cfg.fsample     = 200;
%             cfg.triallength = 500;
%             cfg.ntrials     = 1;
%             cfg.nsignal     = 3;
%             cfg.method      = 'ar';
%             cfg.params = CMat{2,n}; %selected two connections here!
%             cfg.noisecov = NCV;
%             data_raw = ft_connectivitysimulation(cfg);
%             
%             % Do the Signal Mixing
%             data = data_raw;
%             sigmix = repmat(SigMixvec(SMx)/(Nsig-1),Nsig,Nsig).*~eye(Nsig);
%             sigmix = sigmix+eye(Nsig).*1; %(1-NCvec(ncov));
%             data.trial{1} = (data.trial{1} -mean(data.trial{1},2))./std(data.trial{1},[],2);
%             data.trial{1} = sigmix*data.trial{1};
%             
%             % NowDo the ASNR
%             randproc = randn(size(data.trial{1}));
%             for i = 1:size(randproc,1)
%                 s = data.trial{1}(i,:);
%                 s = (s-mean(s))./std(s);
%                 nr  =((SNRvec(i,SNR)*1).*randproc(i,:));
%                 y = s+nr;
%                 snr = var(s)/var(nr);
%                 disp(snr)
% %                 snrbank(SMx,n,i) = snr;
%                 %         y = (y-mean(y))./std(y);
%                 data.trial{1}(i,:) = y;
%             end
%             
%              dataBank{SNR,SMx,n} = data;
%             disp([SNR SMx n])
%         end
%     end
%     mkdir([cd '\benchmark'])
%     save([cd '\benchmark\simdata_10'],'dataBank')
% end

% Now test for recovery with dFC metrics
    load([cd '\benchmark\simdata_10'],'dataBank')

for SNR = 1:size(SNRvec,2)
    for SMx = 1:size(SigMixvec,2)
        bstrap = 1;
        for n = 1:nreps
            
            TrueCMat = CMat{2,n};
            Z = sum(TrueCMat,3);
            Z(Z==0.5) = 0;
            Z(Z==0.3) = 1;
            
            dataN = dataBank{SNR,SMx,n};
            freq = computeSpectra(dataN,[0 0 0],3,0,'-',-1,1);
            [Hz granger grangerft] = computeGranger(freq,cmap(2,:),3,0,'--',1,bstrap);
            [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = computeNPD(dataN,1,8,1,bstrap);
            plotNPD(Hz,npdspctrm,dataN,cmap(3,:),0,':',0)
            if bstrap == 1
                NPG_ci{SNR,SMx}= granger{2,2};
                NPD_ci{SNR,SMx} = npdspctrm{2,2};
                bstrap = 0;
                save([cd '\benchmark\10_NPG_CI_NPD_CI'],'NPG_ci','NPD_ci')
            else
                load([cd '\benchmark\10_NPG_CI_NPD_CI'],'NPG_ci','NPD_ci')
            end
            
            % Now estimate
            NPG = granger{1,2};
            A = squeeze(sum((NPG>NPG_ci{SNR,SMx}),3));
            crit = ceil(size(NPG,3).*0.1);
            Ac = A>crit;
            SA = Ac+Z;
            NPGScore(SNR,SMx,n) = 100.*(sum(SA(:)== 2 | SA(:)== 0) - size(diag(Ac),1))./(numel(Ac) - size(diag(Ac),1));
            
            NPD = npdspctrm{1,2};
            B = squeeze(sum((NPD>NPD_ci{SNR,SMx}),3));
            crit = ceil(size(NPD,3).*0.1);
            Bc = B>crit;
            SB = Bc+Z;
            NPDScore(SNR,SMx,n) = 100.*(sum(SB(:)== 2 | SB(:)== 0) - size(diag(Bc),1))./(numel(Bc) - size(diag(Bc),1));
            %             NPDScore(dataLen,i,n) = sum((B(:)-Z(:)).^2);
            disp([SNR,SMx,n])
        end
    end
end

save([cd '\benchmark\10BenchMarks'],'NPGScore','NPGScore','NCvec','DA')

load([cd '\benchmark\10BenchMarks'],'NPGScore','NPGScore','NCvec','DA')
subplot(1,2,2)
a = plot(NCvec,1-mean(NPGScore,3));
rcmap = brewermap(6,'Reds');
markstylz = {'^','o','s','d'};
for i = 1:ncons
    a(i).Color = rcmap(2+i,:);
    a(i).LineWidth = 2;
    a(i).Marker = markstylz{i};
    a(i).MarkerFaceColor = rcmap(2+i,:);
end
h1 = legend(a,{'1 Connection','2 Connections','3 Connections','4 Connections'},'Location','SouthWest');
grid on
xlabel('Trial Length (2^n)')
ylabel('Estimation Accuracy')
title('Effects of Trial Length on Connection Recovery with mvNPG')
ylim([-3 1.25])
subplot(1,2,1)
b = plot(NCvec,1-mean(NPDScore,3));
bcmap = brewermap(6,'Blues');
for i = 1:ncons
    b(i).Color = bcmap(2+i,:);
    b(i).LineWidth = 2;
    b(i).Marker = markstylz{i};
    b(i).MarkerFaceColor = bcmap(2+i,:);
end
h2 = legend(b,{'1 Connection','2 Connections','3 Connections','4 Connections'},'Location','SouthWest');

grid on
xlabel('Trial Length (2^n)')
ylabel('Estimation Accuracy')
title('Effects of Trial Length on Connection Recovery with NPD')
ylim([-3 1.25])
set(gcf,'Position',[353         512         1120         446])
