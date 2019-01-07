close all; clear all
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools')
%% Figure 8A
load('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\Data\JB\ftdata\cleaned\V6_sources_clean_ROI_OFF_Right_ipsi_SMA.mat')
wrapper_Fig8A_symSNR(vc_clean,3)
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

%% Figure 8B
wrapper_Fig8B_AsymSNR(vc_clean,3)
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

%% Figure 8C
wrapper_Fig8C_sigmix(vc_clean,3)
% figure(2)
titlist = {'STN \leftrightarrow SMA','STN \rightarrow SMA','SMA \rightarrow STN'};
ylimlist = {[-0.05 0.2]; [0 0.1]; [0 0.1]};
for i = 1:3
    subplot(3,1,i); ylim(ylimlist{i});xlim([0 48])
    title(titlist{i})
    legend({'NPD','Granger'})
end
% subplot(N,N,2); ylim([0 0.5])
% subplot(N,N,3); ylim([0 0.5])
set(gcf,'Position',[1263         100         350         979])
