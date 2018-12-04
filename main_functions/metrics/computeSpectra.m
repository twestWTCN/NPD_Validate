function freq = computeSpectra(data,cmap,Nsig,plotfig,linestyle,smord)
if nargin<6
    smord = 2;
end
cfg = [];
cfg.length  =1;
data = ft_redefinetrial(cfg,data);
%% Compute Spectra
cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
cfg.output    = 'fourier';
cfg.tapsmofrq = smord;
freq          = ft_freqanalysis(cfg, data);
if plotfig
    
    for row=1:Nsig
        for col=1:Nsig
            if row==col
                subplot(Nsig,Nsig,(row-1)*Nsig+col);
                Pxy=  mean(abs(squeeze(freq.fourierspctrm(:,col,:))),1);
                [xCalc yCalc b Rsq] = linregress(freq.freq',Pxy');
                Pxy = Pxy-yCalc';
                Pxy = Pxy+abs(min(Pxy));
                plot(freq.freq,Pxy,'color',cmap,'LineWidth',2,'linestyle',linestyle)
                hold on
            end
        end
    end
end