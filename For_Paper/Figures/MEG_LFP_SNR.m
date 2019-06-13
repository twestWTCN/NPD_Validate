close all; clear all
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools')
%% Figure 8A
% load('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\Data\JB\ftdata\cleaned\V6_sources_clean_ROI_OFF_Right_ipsi_SMA.mat')
load('C:\Users\twest\Documents\Work\GitHub\NPD_Validate\data\V6_sources_clean_ROI_ON_Right_ipsi_SMA.mat')

cfg = [];
cfg.bpfilter      = 'yes';
cfg.bpfreq = [20 30];
vc_bp = ft_preprocessing(cfg,vc_clean);

cfg = [];
cfg.bsfilter      = 'yes';
cfg.bsfreq = [20 30];
vc_bs = ft_preprocessing(cfg,vc_clean);


10*log10(var(vc_bp.trial{1},[],2)./var(vc_bs.trial{1},[],2))

