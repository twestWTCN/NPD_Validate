function [Hz granger] = computeGranger(freq,cmap,Nsig,plotfig,linestyle)
cfg           = [];
cfg.method    = 'granger';
grangerft       = ft_connectivityanalysis(cfg, freq);

cfg           = [];
cfg.method    = 'instantaneous_causality';
igrangerft      = ft_connectivityanalysis(cfg, freq);


if plotfig
    for row= 1:Nsig
        for col=1:Nsig
            subplot(Nsig,Nsig,(row-1)*Nsig+col);
            plot(grangerft.freq, squeeze(grangerft.grangerspctrm(row,col,:)),'color',cmap,'LineWidth',2,'linestyle',linestyle)
            hold on
            if col == 1
                ylabel(['signal ' num2str(row)])
            end
            if row == Nsig
                xlabel(['signal ' num2str(col)])
            end
            if row == Nsig && col == floor(Nsig/2)
                xlabel({'to:' ; ['signal ' num2str(col)]})
            end
            if row == floor(Nsig/2) && col == 1
                ylabel({'from:' ; ['signal ' num2str(row)]})
            end
            ylim([0 1])
        end
    end
end


% Put in same format as the NPD
granger{1,1} = igrangerft.instantspctrm;
granger{1,2} = grangerft.grangerspctrm;
granger{1,3} = grangerft.grangerspctrm(2:-1:1,2:-1:1,:);
Hz = grangerft.freq;
