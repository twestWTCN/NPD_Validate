clear all; close all
NC = 4;
% NCvec = linspace(0,1,NC);
NCvec = [0 0.25 0.5 0.75];
for i = 1:NC
    NCtits{i} = ['CV = ' num2str(NCvec(i))];
end

cmap = linspecer(NC);
for ncov = 1:NC
    cfg             = [];
    cfg.ntrials     = 500;
    cfg.triallength = 1;
    cfg.fsample     = 200;
    cfg.nsignal     = 3;
    cfg.method      = 'ar';
    
    cfg.params(:,:,1) = [
        0.5   0   0 ;
        0     0.5   0 ;
        0     0.5     0.5];
    
    cfg.params(:,:,2) = [
       -0.5  0      0 ;
        0   -0.5    0 ;
        0.5  0 -0.5];
    ncv = NCvec(ncov);
    cfg.noisecov      = [ 1    ncv   ncv ;
        ncv    1    ncv ;
        ncv    ncv    1];
    
    data              = ft_connectivitysimulation(cfg);
    
    cfg           = [];
    cfg.method    = 'mtmfft';
    cfg.taper     = 'dpss';
    cfg.output    = 'fourier';
    cfg.tapsmofrq = 2;
    freq          = ft_freqanalysis(cfg, data);
    
    %% COHERENCE
    cfg           = [];
    cfg.method    = 'coh';
    coh           = ft_connectivityanalysis(cfg, freq);
    
    % figure
    % cfg           = [];
    % cfg.parameter = 'cohspctrm';
    % cfg.zlim      = [0 1];
    % ft_connectivityplot(cfg, coh);
    
    figure(1)
    for row=1:3
        for col=1:3
            
            subplot(3,3,(row-1)*3+col);
            plot(coh.freq, squeeze(coh.cohspctrm(row,col,:)),'color',cmap(ncov,:),'LineWidth',2)
            hold on
            if col == 1
                ylabel(['signal ' num2str(row)])
            end
            if row ==3
                xlabel(['signal ' num2str(col)])
            end
            
            ylim([0 1])
        end
    end
    if ncov == NC
        h = legend(NCtits);
        h.Position = [0.7868    0.1352    0.1034    0.1026];
        
        set(gcf,'Position',[680 302.5 836.5 675.5])
        annotation(gcf,'textbox',...
            [0.351 0.941 0.321 0.043],...
            'String',{'Coherence'},...
            'HorizontalAlignment','center',...
            'FitBoxToText','off',...
            'EdgeColor','none');
    end
    %% WPLI
        cfg           = [];
    cfg.method    = 'wpli_debiased';
    wpli           = ft_connectivityanalysis(cfg, freq);
    
    % figure
    % cfg           = [];
    % cfg.parameter = 'cohspctrm';
    % cfg.zlim      = [0 1];
    % ft_connectivityplot(cfg, coh);
    
    figure(2)
    for row=1:3
        for col=1:3
            
            subplot(3,3,(row-1)*3+col);
            plot(wpli.freq, squeeze(wpli.wpli_debiasedspctrm(row,col,:)),'color',cmap(ncov,:),'LineWidth',2)
            hold on
            if col == 1
                ylabel(['signal ' num2str(row)])
            end
            if row ==3
                xlabel(['signal ' num2str(col)])
            end
            
            ylim([0 1])
        end
    end
    if ncov == NC
        h = legend(NCtits);
        h.Position = [0.7868    0.1352    0.1034    0.1026];
        
        set(gcf,'Position',[680 302.5 836.5 675.5])
        annotation(gcf,'textbox',...
            [0.351 0.941 0.321 0.043],...
            'String',{'WPLI'},...
            'HorizontalAlignment','center',...
            'FitBoxToText','off',...
            'EdgeColor','none');
    end
    %% GRANGER
    cfg           = [];
    cfg.method    = 'granger';
    granger       = ft_connectivityanalysis(cfg, freq);
    
    % figure
    % cfg           = [];
    % cfg.parameter = 'grangerspctrm';
    % cfg.zlim      = [0 1];
    % ft_connectivityplot(cfg, granger);
    figure(3)
    for row=1:3
        for col=1:3
            subplot(3,3,(row-1)*3+col);
            plot(granger.freq, squeeze(granger.grangerspctrm(row,col,:)),'color',cmap(ncov,:),'LineWidth',2)
            hold on
            if col == 1
                ylabel(['signal ' num2str(row)])
            end
            if row ==3
                xlabel(['signal ' num2str(col)])
            end
            if row == 2 && col ==1
                ylabel({'from:' ; ['signal ' num2str(row)]})
            end
            ylim([0 1])
        end
    end
    if ncov == NC
        h = legend(NCtits);
        h.Position = [0.7868    0.1352    0.1034    0.1026];
        
        set(gcf,'Position',[680 302.5 836.5 675.5])
        annotation(gcf,'textbox',...
            [0.351 0.941 0.321 0.043],...
            'String',{'Non-Parametric Granger'},...
            'HorizontalAlignment','center',...
            'FitBoxToText','off',...
            'EdgeColor','none');
    end
    % cfg = [];
    % cfg.viewmode = 'butterfly';  % you can also specify 'butterfly'
    % ft_databrowser(cfg, data);
    
    for i = 1:length(data.label)
        for j = 1:length(data.label)
            %                 ftspect = freqFour.fourierspctrm(:,[chcombs(1,i) chcombs(2,i)],:);
            %                 [f13,t13,cl13] = sp2a2_R2_TW(FTdata.fsample,8,ftspect);
            x = []; y = [];
            for p = 1:size(data.trial,2)
                x = [x data.trial{p}(i,:)];
                y = [y data.trial{p}(j,:)];
            end
            %                 [f13,~,~]=sp2a2_R2_mt(x',y',FTdata.fsample,7,'M1');
            [f13,t13,c13]=sp2a2_R2(x',y',data.fsample,8);
            npdspctrm{1,1}(i,j,:) = f13(:,10);
            npdspctrm{1,2}(i,j,:) = f13(:,12); % Backward
            npdspctrm{1,3}(i,j,:) = f13(:,11); % Forward
            
            npdcrcv(i,j,:) = t13(:,3);
            nscohspctrm(i,j,:) = f13(:,4);
        end
    end
    % Partialised
    for i = 1:length(data.label)
        for j = 1:length(data.label)
            %                 ftspect = freqFour.fourierspctrm(:,[chcombs(1,i) chcombs(2,i)],:);
            %                 [f13,t13,cl13] = sp2a2_R2_TW(FTdata.fsample,8,ftspect);
            x = []; y = []; z = [];
            for p = 1:size(data.trial,2)
                x = [x data.trial{p}(i,:)];
                y = [y data.trial{p}(j,:)];
                z = [z data.trial{p}(1,:)];
            end
            %                 [f13,~,~]=sp2a2_R2_mt(x',y',FTdata.fsample,7,'M1');
            [f13,t13,c13]=sp2_R2a_pc1(x,y,z,data.fsample,2^8);
            npdspctrmZ{1,1}(i,j,:) = f13(:,10);
            npdspctrmZ{1,2}(i,j,:) = f13(:,12); % Backward
            npdspctrmZ{1,3}(i,j,:) = f13(:,11); % Forward
            
            npdcrcvZ(i,j,:) = t13(:,3);
            nscohspctrmZ(i,j,:) = f13(:,4);
        end
    end
    %% PLOT NPD
    figure(4)
    for i = 1:length(data.label)
        for j = 1:length(data.label)
            linind = sub2ind([length(data.label),length(data.label)],j,i);
            subplot(length(data.label),length(data.label),linind)
            if i>j
                fy = squeeze(npdspctrm{1,2}(j,i,:));
            elseif j>i
                fy = squeeze(npdspctrm{1,2}(j,i,:));
            elseif j==i
                fy = squeeze(npdspctrm{1,1}(i,j,:));
            end
            fx = f13(:,1);
            
            plot(fx,fy,'color',cmap(ncov,:),'LineWidth',2)
            hold on
            if j == 1
                ylabel(['signal ' num2str(i)])
            end
            if i ==3
                xlabel(['signal ' num2str(j)])
            end
            if i == 2 && j ==1
                ylabel({'from:' ; ['signal ' num2str(i)]})
            end
            if ncov == NC
                plot(fx,repmat(c13.ch_c95,1,size(fy,1)),'k--')
            end
            ylim([0 1])
            clear fy
        end
    end
    if ncov==NC
        h = legend(NCtits);
        h.Position = [0.7868    0.1352    0.1034    0.1026];
        
        set(gcf,'Position',[680 302.5 836.5 675.5])
        annotation(gcf,'textbox',...
            [0.351 0.941 0.321 0.043],...
            'String',{'NonParametric Directionality'},...
            'HorizontalAlignment','center',...
            'FitBoxToText','off',...
            'EdgeColor','none');
    end
%% NPD CROSS CORRELATION    
        figure(5)
    for i = 1:length(data.label)
        for j = 1:length(data.label)
            linind = sub2ind([length(data.label),length(data.label)],j,i);
            subplot(length(data.label),length(data.label),linind)
            fy = squeeze(npdcrcv(i,j,:));
            fx = t13(:,1);
            
            plot(fx,fy,'color',cmap(ncov,:),'LineWidth',2)
            hold on
            if j == 1
                ylabel(['signal ' num2str(i)])
            end
            if i ==3
                xlabel(['signal ' num2str(j)])
            end
            if i == 2 && j ==1
                ylabel({'from:' ; ['signal ' num2str(i)]})
            end
            if ncov == NC
                plot(fx,repmat(c13.rho_c95,1,size(fy,1)),'k--')
                hold on
                plot(fx,repmat(-c13.rho_c95,1,size(fy,1)),'k--')
            end
            ylim([-0.5 1]);
            xlim([-50 50])
            clear fy
        end
    end
    if ncov==NC
        h = legend(NCtits);
        h.Position = [0.7868    0.1352    0.1034    0.1026];
        
        set(gcf,'Position',[680 302.5 836.5 675.5])
        annotation(gcf,'textbox',...
            [0.351 0.941 0.321 0.043],...
            'String',{'NonParametric Directionality RhoYX'},...
            'HorizontalAlignment','center',...
            'FitBoxToText','off',...
            'EdgeColor','none');
    end
    %% PLOT NPDZ
    figure(7)
    for i = 1:length(data.label)
        for j = 1:length(data.label)
            linind = sub2ind([length(data.label),length(data.label)],j,i);
            subplot(length(data.label),length(data.label),linind)
            if i>j
                fy = squeeze(npdspctrmZ{1,2}(j,i,:));
            elseif j>i
                fy = squeeze(npdspctrmZ{1,2}(j,i,:));
            elseif j==i
                fy = squeeze(npdspctrmZ{1,1}(i,j,:));
            end
            fx = f13(:,1);
            
            plot(fx,fy,'color',cmap(ncov,:),'LineWidth',2)
            hold on
            if j == 1
                ylabel(['signal ' num2str(i)])
            end
            if i ==3
                xlabel(['signal ' num2str(j)])
            end
            if i == 2 && j ==1
                ylabel({'from:' ; ['signal ' num2str(i)]})
            end
            if ncov == NC
                plot(fx,repmat(c13.ch_c95,1,size(fy,1)),'k--')
            end
            ylim([0 1])
            clear fy
        end
    end
    if ncov==NC
        h = legend(NCtits);
        h.Position = [0.7868    0.1352    0.1034    0.1026];
        
        set(gcf,'Position',[680 302.5 836.5 675.5])
        annotation(gcf,'textbox',...
            [0.351 0.941 0.321 0.043],...
            'String',{'NonParametric Directionality Partialed on 3'},...
            'HorizontalAlignment','center',...
            'FitBoxToText','off',...
            'EdgeColor','none');
    end
%% NPDZ CROSS CORRELATION    
        figure(8)
    for i = 1:length(data.label)
        for j = 1:length(data.label)
            linind = sub2ind([length(data.label),length(data.label)],j,i);
            subplot(length(data.label),length(data.label),linind)
            fy = squeeze(npdcrcvZ(i,j,:));
            fx = t13(:,1);
            
            plot(fx,fy,'color',cmap(ncov,:),'LineWidth',2)
            hold on
            if j == 1
                ylabel(['signal ' num2str(i)])
            end
            if i ==3
                xlabel(['signal ' num2str(j)])
            end
            if i == 2 && j ==1
                ylabel({'from:' ; ['signal ' num2str(i)]})
            end
            if ncov == NC
                plot(fx,repmat(c13.rho_c95,1,size(fy,1)),'k--')
                hold on
                plot(fx,repmat(-c13.rho_c95,1,size(fy,1)),'k--')
            end
            ylim([-0.5 1]);
            xlim([-50 50])
            clear fy
        end
    end
    if ncov==NC
        h = legend(NCtits);
        h.Position = [0.7868    0.1352    0.1034    0.1026];
        
        set(gcf,'Position',[680 302.5 836.5 675.5])
        annotation(gcf,'textbox',...
            [0.351 0.941 0.321 0.043],...
            'String',{'NonParametric Directionality RhoYX Partialed'},...
            'HorizontalAlignment','center',...
            'FitBoxToText','off',...
            'EdgeColor','none');
    end

end
saveallfiguresSFLAP('C:\Users\Tim\Documents\Work\LitvakProject\NPD_Verify\MVAR_model\figures\1link\outputs_1link','-jpg')