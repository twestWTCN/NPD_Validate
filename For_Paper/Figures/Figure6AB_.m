% NPD_Validate_AddPaths()
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools')
clear; close all
for X = 1:2
    fresh = 0; % (1) Run permutation tests; (0) Load precomputed tables.
    permrun = X; % (1) FFT Shuffle; (2) Phase Randomize
    %% This script will reproduce Figure 6A and B -
    % Investigating the role of data availability upon the accuracy
    % of connectivity recovery when using non-parametric directionality
    % (NPD), and non-parametric Granger causality (NPG).
    
    % Plot Colors
    cmap = linspecer(4);
    plotfig = 0;
    % Setup random Graphs
    rng(231284)
    ncons = 3;
    nreps = 24;
    [CMat,NCV] = makeRndGraphs(ncons,nreps,3);
    NCvec = linspace(3,10,8); %linspace(3,12,11);
    DA = 500; % Total data availability (200s)
    fsamp= 200;
    Nsig = 3;
    %% Simulate data
        if fresh == 1
            for dataLen = 1:size(NCvec,2)
                clear data
                rng(63487)
                for i = 1:ncons
                    parfor n = 1:nreps
                        %             Simulate Data
                        cfg             = [];
                        cfg.fsample     = fsamp;
                        cfg.triallength = (2.^(NCvec(dataLen)))./cfg.fsample;
                        cfg.ntrials     = floor(DA./cfg.triallength);
                        cfg.nsignal     = Nsig;
                        cfg.method      = 'ar';
                        cfg.params = CMat{i,n};
                        cfg.noisecov = NCV;
                        data{i,n} = ft_connectivitysimulation(cfg);
                        disp([dataLen i n])
                    end
                end
                mkdir([cd '\benchmark'])
                save([cd '\benchmark\simdata_9A_' num2str(dataLen)],'data','CMat')
            end
        end
    NPG_ci = {}; NPG_ci = {};
    
    %% Now test for recovery with dFC metrics
    for dataLen = 1:size(NCvec,2)
        load([cd '\benchmark\simdata_9A_' num2str(dataLen)],'data','CMat')
        
        if fresh == 1
            perm = 1; permtype = permrun;
        else
            perm = 0; permtype = permrun;
            load([cd '\benchmark\9A_NPG_CI_NPD_CI_' num2str(permrun)],'NPG_ci','NPD_ci')
        end
        
        for i = 1:ncons
            dataN = data(i,:);
            parfor n = 1:nreps
                TrueCMat = CMat{i,n};
                Z = sum(TrueCMat,3);
                Z(Z==0.5) = 0; % These are the self connections
                Z = Z>0;
                dataIN = dataN{n};
                % Compute spectra
                datalength =(2.^(NCvec(dataLen)))./fsamp;
                freq = computeSpectra(dataIN,[0 0 0],Nsig,plotfig,'-',-1,[]);
                % Compute connectivity
                [Hz granger grangerft] = computeGranger(freq,1,perm,permtype)
                [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = ft_computeNPD(freq,fsamp,1,NCvec(dataLen),perm,permtype);
%                                  plotNPD(Hz,npdspctrm,dataIN,cmap(3,:),1,'-',0)
%                                  plotNPD(grangerft.freq,granger,freq,cmap(2,:),1,'-',0)
%                 if perm == 1
%                     NPG_ci{dataLen}= granger{2,2};
%                     NPD_ci{dataLen} = npdspctrm{2,2};
%                     perm = 0;
%                     save([cd '\benchmark\9A_NPG_CI_NPD_CI_' num2str(permrun)],'NPG_ci','NPD_ci')
%                 else
%                     load([cd '\benchmark\9A_NPG_CI_NPD_CI_' num2str(permrun)],'NPG_ci','NPD_ci')
%                 end
                
                % Now estimate
                NPG = granger{1,2};
                CT = nanmean(NPG_ci{12}(:));
                A = squeeze(sum((NPG>CT),3));
                crit = ceil(size(NPG,3).*0.1);
                Ac = A>crit;
                NPGScore(dataLen,i,n) = matrixScore(Ac,Z);
                
                NPD = npdspctrm{1,2};
                CT = mean(NPD_ci{12}(:));
                B = squeeze(sum((NPD>CT),3));
                crit = ceil(size(NPD,3).*0.1);
                Bc = B>crit;
                NPDScore(dataLen,i,n) = matrixScore(Bc,Z);
                %             NPDScore(dataLen,i,n) = sum((B(:)-Z(:)).^2);
                disp([n i dataLen])
            end
        end
    end
    % 2
    save([cd '\benchmark\9ABenchMarks_' num2str(permrun)],'NPDScore','NPGScore','NCvec','DA')
end

for X = 1:2
    load([cd '\benchmark\9ABenchMarks_' num2str(X)],'NPDScore','NPGScore','NCvec','DA')
    
    figure(X)
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
%     h1 = legend(a,{'1 Connection','2 Connections','3 Connections','4 Connections'},'Location','SouthWest');
    grid on
    xlabel('Trial Length ') 
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
%     h2 = legend(b,{'1 Connection','2 Connections','3 Connections','4 Connections'},'Location','SouthWest');
    
    grid on
    xlabel('Trial Length ')
    ylabel('Estimation Accuracy')
    title('Effects of Trial Length on Connection Recovery with NPD')
    set(gcf,'Position',[353         512         1120         446])
    ylim([0 100])
end