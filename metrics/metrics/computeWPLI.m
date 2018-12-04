function computeWPLI(freq,cmap,Nsig,plotfig,linespec)

    cfg           = [];
    cfg.method    = 'wpli_debiased';
    wpli           = ft_connectivityanalysis(cfg, freq);
    
    % figure
    % cfg           = [];
    % cfg.parameter = 'cohspctrm';
    % cfg.zlim      = [0 1];
    % ft_connectivityplot(cfg, coh);
    
    figure(6)
    for row=1:Nsig
        for col=1:Nsig
            
            subplot(Nsig,Nsig,(row-1)*Nsig+col);
            plot(wpli.freq, squeeze(wpli.wpli_debiasedspctrm(row,col,:)),'color',cmap,'LineWidth',2,'linestyle',linestyle)
            hold on
            if col == 1
                ylabel(['signal ' num2str(row)])
            end
            if row ==Nsig
                xlabel(['signal ' num2str(col)])
            end
            
            ylim([0 1])
        end
    end
    if plotfig
        set(gcf,'Position',[636 172 1264 889])
        h = legend(NCtits);
        set(h,'Position',[0.895 0.845 0.0837 0.0804],'FontSize',8)
        annotation(gcf,'textbox',...
            [0.351 0.941 0.321 0.043],...
            'String',{'WPLI'},...
            'HorizontalAlignment','center',...
            'FitBoxToText','off',...
            'EdgeColor','none');
    end
