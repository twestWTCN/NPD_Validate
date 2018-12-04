function plotNSCoherence(coh,cmap,Nsig,plotfig,linestyle)
if plotfig
    for row=1:Nsig
        for col=1:Nsig
            
            subplot(Nsig,Nsig,(row-1)*Nsig+col);
            plot(coh.freq, squeeze(coh.cohspctrm(row,col,:)),'color',cmap,'LineWidth',2,'linestyle',linestyle)
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
end