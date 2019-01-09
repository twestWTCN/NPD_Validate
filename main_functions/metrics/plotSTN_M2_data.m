function plotSTN_M2_data(Hz,npdspctrm,data,cmap,plotfig,linestyle,confindflag)
%% PLOT NPD
if plotfig
    for i = 1:3
        subplot(3,1,i)
        if i ==1
            fy = squeeze(npdspctrm{1,1}(2,1,:)); % i = 1 instant;
            cy = squeeze(npdspctrm{2,2}(2,1,:));
        elseif i == 2
            fy = squeeze(npdspctrm{1,2}(1,2,:)); % i = 2 reverse (j->i SMA->STN);
            cy = squeeze(npdspctrm{2,2}(1,2,:));
        elseif i ==3
            fy = squeeze(npdspctrm{1,2}(2,1,:)); % i = 3 forward (i->j STN -> SMA);
            cy = squeeze(npdspctrm{2,1}(2,1,:));
        end
        fx = Hz;
        plot(fx,fy,'color',cmap,'LineWidth',2,'linestyle',linestyle)
        hold on
        if confindflag == 1
            plot(fx,cy,'color',cmap,'linestyle',':','LineWidth',1.5)
        end
        
        %         if plotfig
        %             plot(fx,repmat(c13.ch_c95,1,size(fy,1)),'k--')
        %         end
        ylim([0 0.1]); grid on
        clear fy
    end
end

