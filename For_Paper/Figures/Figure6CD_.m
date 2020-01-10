% NPD_Validate_AddPaths()
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools')
close all

%% This script will reproduce Figure 6C and D (-
% Investigating the role of data availability upon the accuracy
% of connectivity recovery when using non-parametric directionality
% (NPD), and non-parametric Granger causality (NPG).

% Plot Colors
cmap = linspecer(4);

% Setup random Graphs
rng(42342)
ncons = 4;
nreps = 20; %25
sz = 3;
[CMat,NCV] = makeRndGraphs(ncons,nreps,sz);
NCvec = linspace(3,12,12); %12
DT = 35;
fsamp= 200;
Nsig = 3;

%% Simulate data
for dataLen = 1:size(NCvec,2)
    clear data
    rng(63487)
    for i = 1:ncons
        for n = 1:nreps
            %             Simulate Data
            cfg             = [];
            cfg.fsample     = fsamp;
            cfg.triallength = (2.^(NCvec(dataLen)))./cfg.fsample;
            cfg.ntrials     = DT;
            cfg.nsignal     = Nsig;
            cfg.method      = 'ar';
            cfg.params = CMat{i,n};
            cfg.noisecov = NCV;
            data{i,n} = ft_connectivitysimulation(cfg);
            disp([dataLen i n])
        end
    end
    mkdir([cd '\benchmark'])
    save([cd '\benchmark\simdata_9B_' num2str(dataLen)],'data')
end

%% Now test for recovery with dFC metrics
for dataLen = 1:numel(NCvec)
    load([cd '\benchmark\simdata_9B_' num2str(dataLen)],'data')
    perm = 1; permtype = 2;
    for i = 1:ncons
        load([cd '\benchmark\9B_NPG_CI_NPD_CI'],'NPG_ci','NPD_ci')
        for n = 1:nreps
            TrueCMat = CMat{i,n};
            Z = sum(TrueCMat,3);
            Z(Z==0.5) = 0;
            Z(Z==0.3) = 1;
            
            dataN = data{i,n};
            datalength =(2.^(NCvec(dataLen)))./fsamp;
            freq = computeSpectra(dataN,[0 0 0],3,0,'-',-1,dataN.cfg.triallength);
            [Hz granger grangerft] = computeGranger(freq,cmap(2,:),3,1,'--',1,perm);
            [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = ft_computeNPD(freq,fsamp,1,NCvec(dataLen),perm,permtype);
            
            if perm== 1
                NPG_ci{dataLen}= granger{2,2};
                NPD_ci{dataLen} = npdspctrm{2,2};
                perm = 0;
                save([cd '\benchmark\9B_NPG_CI_NPD_CI'],'NPG_ci','NPD_ci')
            else
                load([cd '\benchmark\9B_NPG_CI_NPD_CI'],'NPG_ci','NPD_ci')
            end
            
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
            NPDScore(dataLen,i,n) = matrixScore(Bc,Z);
            disp([n i dataLen])
            
        end
    end
end

save([cd '\benchmark\9BBenchMarks'],'NPDScore','NPGScore','NCvec')

load([cd '\benchmark\9BBenchMarks.mat'],'NPDScore','NPGScore','NCvec')
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
% ylim([-3 1.25])
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
% ylim([-3 1.25])
set(gcf,'Position',[353         512         1120         446])
