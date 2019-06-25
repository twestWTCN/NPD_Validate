close all; clear
addpath('C:\Users\Tim\Documents\Work\GIT\BrewerMap')
cmap = linspecer(4);

%% Simulation 9A - Data Length with Fixed Data Length (25s)
rng(927693)
ncons = 4;
nreps = 20;
[CMat,NCV] = makeRndGraphs(ncons,nreps,3);
NCvec = linspace(5,12,12);
DA = 100; % Total data availability (100s)

% Nsamps = 150;

% Simulate data
for dataLen = 1:size(NCvec,2)
    clear data
    rng(7865734)
    for i = 1:ncons
        parfor n = 1:nreps
%             Simulate Data
            cfg             = [];
            cfg.fsample     = 200;
            cfg.triallength = (2.^(NCvec(dataLen)))./cfg.fsample;
            cfg.ntrials     = ceil(DA./cfg.triallength);
            cfg.nsignal     = 3;
            cfg.method      = 'ar';
            cfg.params = CMat{i,n};
            cfg.noisecov = NCV;
            data{i,n} = ft_connectivitysimulation(cfg);
            disp([dataLen i n])
        end
    end
    mkdir([cd '\benchmark'])
    save([cd '\benchmark\simdata_9A_' num2str(dataLen)],'data')
end

% Now test for recovery with dFC metrics
for dataLen = 1:size(NCvec,2)
    load([cd '\benchmark\simdata_9A_' num2str(dataLen)],'data')
    bstrap = 0;
    for i = 1:ncons
        load([cd '\benchmark\9A_NPG_CI_NPD_CI'],'NPG_ci','NPD_ci')
        for n = 1:nreps
            TrueCMat = CMat{i,n};
            Z = sum(TrueCMat,3);
            Z(Z==0.5) = 0;
            Z(Z==0.3) = 1;
            
            dataN = data{i,n};
            freq = computeSpectra(dataN,[0 0 0],3,0,'-',-1,dataN.cfg.triallength);
            [Hz granger grangerft] = computeGranger(freq,cmap(2,:),3,0,'--',1,bstrap);
            [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = computeNPD(dataN,1,(NCvec(dataLen)),1,bstrap);
            %             plotNPD(Hz,npdspctrm,dataN,cmap(3,:),0,':',0)
%             if bstrap == 1
%                 NPG_ci{dataLen}= granger{2,2};
%                 NPD_ci{dataLen} = npdspctrm{2,2};
%                 bstrap = 0;
%                 save([cd '\benchmark\9A_NPG_CI_NPD_CI'],'NPG_ci','NPD_ci')
%             else
%                 load([cd '\benchmark\9A_NPG_CI_NPD_CI'],'NPG_ci','NPD_ci')
%             end
            
            % Now estimate
            NPG = granger{1,2};
            A = squeeze(sum((NPG>NPG_ci{dataLen}),3));
            crit = ceil(size(NPG,3).*0.1);
            Ac = A>crit;
            NPGScore(dataLen,i,n) = matrixScore(Ac,Z);
            
            NPD = npdspctrm{1,2};
            B = squeeze(sum((NPD>NPD_ci{dataLen}),3));
            crit = ceil(size(NPD,3).*0.1);
            Bc = B>crit;
            NPDScore(dataLen,i,n) = matrixScore(Bc,Z);1--
            %             NPDScore(dataLen,i,n) = sum((B(:)-Z(:)).^2);
            disp([n i dataLen])
        end
    end
end
% 2
save([cd '\benchmark\9ABenchMarks'],'NPDScore','NPGScore','NCvec','DA')

load([cd '\benchmark\9ABenchMarks'],'NPDScore','NPGScore','NCvec','DA')
subplot(1,2,2)
a = plot(NCvec,mean(NPGScore,3));
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
ylim([0 100])
subplot(1,2,1)
b = plot(NCvec,mean(NPDScore,3));
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
set(gcf,'Position',[353         512         1120         446])
ylim([0 100])
