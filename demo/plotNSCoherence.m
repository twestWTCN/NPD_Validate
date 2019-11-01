function plotNSCoherence(coh,cmap,Nsig,plotfig,linestyle,confindflag)
if nargin<6
    confindflag = 0;
end
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
        if confindflag == 1 && i~=j
            plot(coh.freq,squeeze(coh.ci(row,col,:)),'color',cmap,'linestyle',':','LineWidth',1.5)
        end
        ylim([0 1])
        end
    end
end