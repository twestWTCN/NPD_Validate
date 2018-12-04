function plotNPD_zero(Hz,npdspctrm,data,cmap,plotfig,linestyle)
%% PLOT NPD
if plotfig
for i = 1:length(data.label)
    for j = 1:length(data.label)
        linind = sub2ind([length(data.label),length(data.label)],j,i);
        subplot(length(data.label),length(data.label),linind)
        
        if i>j
            fy = squeeze(npdspctrm{1,1}(j,i,:));
        elseif j>i
            fy = squeeze(npdspctrm{1,1}(j,i,:));
        elseif j==i
            fy = squeeze(npdspctrm{1,1}(i,j,:));
        end
        fx = Hz;
        
        plot(fx,fy,'color',cmap,'LineWidth',2,'linestyle',linestyle)
        hold on
        if j == 1
            ylabel(['signal ' num2str(i)])
        end
        if i == length(data.label)
            xlabel(['signal ' num2str(j)])
        end
        if i == 3 && j == 2
            xlabel({'to::' ; ['signal ' num2str(i)]})
        end
        if i == 2 && j == 1
            ylabel({'from:' ; ['signal ' num2str(i)]})
        end
%         if plotfig
%             plot(fx,repmat(c13.ch_c95,1,size(fy,1)),'k--')
%         end
        ylim([0 1]); grid on
        clear fy
    end
end
end
