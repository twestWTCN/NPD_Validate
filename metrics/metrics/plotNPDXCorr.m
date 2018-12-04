function plotNPDXCorr(lags,npdcrcvZ,data,cmap,plotfig,linestyle)
    for i = 1:length(data.label)
        for j = 1:length(data.label)
            linind = sub2ind([length(data.label),length(data.label)],j,i);
            subplot(length(data.label),length(data.label),linind)
            fy = squeeze(npdcrcvZ(i,j,:));
            fx = lags;
            
            plot(fx,fy,'color',cmap,'LineWidth',2,'linestyle',linestyle)
            hold on
            if j == 1
                ylabel({['signal ' num2str(i)];'R(x)'})
            end
            if i ==length(data.label)
                xlabel({'lag (s)';['signal ' num2str(j)]})
            end
            if i == 2 && j ==1
                ylabel({'from:' ; ['signal ' num2str(i)];'R(x)'})
            end
            
             if i == 3 && j ==2
                xlabel({'lag (s)';['signal ' num2str(j)]; 'to:'})
            end
            if plotfig
                plot(fx,repmat(c13.rho_c95,1,size(fy,1)),'k--')
                hold on
                plot(fx,repmat(-c13.rho_c95,1,size(fy,1)),'k--')
            end
            ylim([-0.5 1]);
            xlim([-500 500]); grid on
            clear fy
        end
    end
    if plotfig
        set(gcf,'Position',[636 172 1264 889])
        h = legend(NCtits);
        set(h,'Position',[0.895 0.845 0.0837 0.0804],'FontSize',8)
        annotation(gcf,'textbox',...
            [0.351 0.941 0.321 0.043],...
            'String',{'NonParametric Directionality RhoYX Partialed on ' num2str(frstord) ': XCOR'},...
            'HorizontalAlignment','center',...
            'FitBoxToText','off',...
            'EdgeColor','none');
    end
