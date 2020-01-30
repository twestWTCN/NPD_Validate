% NPD_Validate_AddPaths()
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools')
close all
fresh = 1; % (1) Run permutation tests; (0) Load precomputed tables.
permrun = 2; % (1) FFT Shuffle; (2) Phase Randomize

%% This script will reproduce Figure 9 - 
% Testing for the confounding effects of symmetric and asymmetric SNR,
% and instantaneous signal mixing upon estimation of directed functional
% connectivity in experimental data recorded in patients with Parkinson’s disease. 

% Load in the data

% load('C:\Users\twest\Documents\Work\GitHub\NPD_Validate\data\V6_sources_clean_ROI_OFF_Right_ipsi_SMA.mat'); % FREE
load('C:\Users\Tim\Documents\Work\GIT\NPD_Validate\data\V6_sources_clean_ROI_OFF_Right_ipsi_SMA.mat'); % SFLAP
% load('C:\Users\timot\Documents\GitHub\NPD_Validate\data\V6_sources_clean_ROI_OFF_Right_ipsi_SMA.mat'); % SFLAP

%% Figure 9A
% Modify SNR Symmetrically
figure(1)
wrapper_Fig9A_symSNR(vc_clean,3,permrun,fresh)

titlist = {'STN \leftrightarrow SMA','STN \rightarrow SMA','SMA \rightarrow STN'};
ylimlist = {[0 0.15]; [0 0.15]; [0 0.15]};
for i = 1:3
    subplot(3,1,i); ylim(ylimlist{i});xlim([0 48])
    title(titlist{i})
    legend({'NPD','Granger'})
end
% subplot(N,N,2); ylim([0 0.5])
% subplot(N,N,3); ylim([0 0.5])
set(gcf,'Position',[1263         100         350         979])

%% Figure 9B
% Modify SNR Asymmetrically
figure(2)
wrapper_Fig9B_AsymSNR(vc_clean,3,permrun,fresh)

titlist = {'STN \leftrightarrow SMA','STN \rightarrow SMA','SMA \rightarrow STN'};
ylimlist = {[0 0.15]; [0 0.15]; [0 0.15]};
for i = 1:3
    subplot(3,1,i); ylim(ylimlist{i});xlim([0 48])
    title(titlist{i})
    legend({'NPD','Granger'})
end
% subplot(N,N,2); ylim([0 0.5])
% subplot(N,N,3); ylim([0 0.5])
set(gcf,'Position',[1263         100         350         979])

%% Figure 8C
% Modify Signal Mixing
figure(3)
wrapper_Fig9C_sigmix(vc_clean,3,permrun,fresh)

titlist = {'STN \leftrightarrow SMA','STN \rightarrow SMA','SMA \rightarrow STN'};
ylimlist = {[-0.1 0.2]; [0 0.15]; [0 0.15]};
for i = 1:3
    subplot(3,1,i); ylim(ylimlist{i});xlim([0 48])
    title(titlist{i})
    legend({'NPD','Granger'})
end
% subplot(N,N,2); ylim([0 0.5])
% subplot(N,N,3); ylim([0 0.5])
set(gcf,'Position',[1263         100         350         979])
