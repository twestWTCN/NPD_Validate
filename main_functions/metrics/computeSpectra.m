function freq = computeSpectra(data,cmap,Nsig,plotfig,linestyle,smord,dl)
if nargin<6
    smord = -1;
end
if nargin<7
    dl = 1;
end
if smord == -1
    tname = 'hanning';
else
    tname = 'dpss';
end
if ~isempty(dl)
    cfg = [];
    cfg.length  = dl;
    data = ft_redefinetrial(cfg,data);
end
%% Compute Spectra
cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = tname;
cfg.output    = 'fourier';
if smord >0
    cfg.tapsmofrq = smord;
end
cfg.keeptrials = 'yes';
cfg.padtype = 'zero';
% cfg.pad = 2;
freq          = ft_freqanalysis(cfg, data);

if plotfig
    for row=1:Nsig
        for col=1:Nsig
            if row==col
                subplot(Nsig,Nsig,(row-1)*Nsig+col);
                Pxy=  mean(abs(squeeze(freq.fourierspctrm(:,col,:))),1);
                %                 [xCalc yCalc b Rsq] = linregress(freq.freq',Pxy');
                %                 Pxy2 = Pxy-yCalc';
                %                 Pxy2 = Pxy+abs(min(Pxy2));
                plot(freq.freq,Pxy,'color',cmap,'LineWidth',2,'linestyle',linestyle)
                hold on
            end
        end
    end
end