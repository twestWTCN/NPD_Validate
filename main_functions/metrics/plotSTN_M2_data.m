function plotSTN_M2_data(Hz,npdspctrm,data,cmap,plotfig,linestyle)
%% PLOT NPD
if plotfig
    for i = 1:3
        subplot(3,1,i)
        fy = squeeze(npdspctrm{1,i}(1,2,:));
        fx = Hz;
        plot(fx,fy,'color',cmap,'LineWidth',2,'linestyle',linestyle)
        hold on
        
        
        %         if plotfig
        %             plot(fx,repmat(c13.ch_c95,1,size(fy,1)),'k--')
        %         end
        ylim([0 0.1]); grid on
        clear fy
    end
end

