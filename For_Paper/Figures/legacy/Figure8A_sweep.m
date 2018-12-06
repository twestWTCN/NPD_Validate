close all; clear all
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools')
%% Simulation 8- STN/M2 Data
R.subname = {'DF','DP','DS','JB','JN','LM','LN02','LN03','MC','WB','OXSH_D12'}; % JA,,LN01: Bad Recordings - MW borderline - Rhs is ok
R.condname = {'OFF','ON'};
R.siden = {'Left','Right'};

for sub = 4; %3:length(R.subname)
    for cond = 1; %:length(R.condname)
        for side = 2
            load(['C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\Data\' R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_ipsi_SMA.mat'])
            mvarconsim_npdver_F8A(vc_clean,1)
            % figure(2)
            titlist = {'STN \leftrightarrow SMA','STN \rightarrow SMA','SMA \rightarrow STN'};
            ylimlist = {[0 0.1]; [0 0.1]; [0 0.1]};
            for i = 1:3
                subplot(3,1,i); ylim(ylimlist{i});xlim([0 48])
                title(titlist{i})
                legend({'NPD','Granger'})
            end
            % subplot(N,N,2); ylim([0 0.5])
            % subplot(N,N,3); ylim([0 0.5])
            set(gcf,'Position',[1263         100         350         979])
            [sub side cond]
            pause(5)
            close all

        end
    end
end



