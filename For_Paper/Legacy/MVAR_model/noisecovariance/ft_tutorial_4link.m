clear all; close all
NC = 5;
% NCvec = linspace(0,1,NC);
NCvec = [0 0.1 0.3 0.5 0.7];
for i = 1:NC
    NCtits{i} = ['SigLeak = ' num2str(NCvec(i))];
end

cmap = linspecer(NC);
for ncov = 1:NC
    ncv = NCvec(ncov);
    Nsig = 4;
    cfg             = [];
    cfg.ntrials     = 1;
    cfg.triallength = 500;
    cfg.fsample     = 200;
    cfg.nsignal     = Nsig;
    cfg.method      = 'ar';
    
    cfg.params(:,:,1) = [
        0.5   0     0    0   ;
        0   0.5     0    0   ;
        0     0     0.5  0   ;
        0     0     0    0.5 ];
    cfg.params(2,1,1) = 0.5;
    cfg.params(3,1,1) = 0.5;
    cfg.params(4,2,1) = 0.5;
    
    cfg.params(:,:,2) = [
        -0.5   0     0    0   ;
        0    -0.5   0    0   ;
        0     0   -0.5   0   ;
        0     0     0   -0.5 ];
    
%     cfg.params(:,:,3) = [
%         0.5   0     0    0   ;
%         0     0.5   0    0   ;
%         0     0     0.5  0   ;
%         0     0     0   0.5 ];
    
    
    cfg.noisecov      = eye(Nsig);
    
    data              = ft_connectivitysimulation(cfg);
    
    % Weighted mixture of signals
    %     sigmix = [
    %         0.7 0.1 0.1 0.1;
    %         0.1 0.7 0.1 0.1;
    %         0.1 0.1 0.7 0.1;
    %         0.1 0.1 0.1 0.7;
    %         ];
    sigmix = repmat(NCvec(ncov)/(Nsig-1),Nsig,Nsig).*~eye(Nsig);
    sigmix = sigmix+eye(Nsig).*(1-NCvec(ncov));
    
    data.trial{1} = sigmix*data.trial{1};
    
    %     sigmix = [1 1 1 1];
    %     volcon = [.1 .1 .1 .1];
    %     for i = 1:NC
    %         x(i,:) = data.trial{1}(i,:).*sigmix(i);
    %     end
    %     x = mean(x,1);
    %     for i = 1:NC
    %         data.trial{1}(i,:) = data.trial{1}(i,:) + volcon(i).*x;
    %     end
    cfg = [];
    cfg.length  =1;
    data = ft_redefinetrial(cfg,data);
    
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
    for row=1:Nsig
        for col=1:Nsig
            
            subplot(Nsig,Nsig,(row-1)*Nsig+col);
            plot(coh.freq, squeeze(coh.cohspctrm(row,col,:)),'color',cmap(ncov,:),'LineWidth',2)
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
    if ncov == NC
        set(gcf,'Position',[636 172 1264 889])
        h = legend(NCtits);
        set(h,'Position',[0.895 0.845 0.0837 0.0804],'FontSize',8)
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
    for row=1:Nsig
        for col=1:Nsig
            
            subplot(Nsig,Nsig,(row-1)*Nsig+col);
            plot(wpli.freq, squeeze(wpli.wpli_debiasedspctrm(row,col,:)),'color',cmap(ncov,:),'LineWidth',2)
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
    if ncov == NC
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
    for row=1:Nsig
        for col=1:Nsig
            subplot(Nsig,Nsig,(row-1)*Nsig+col);
            plot(granger.freq, squeeze(granger.grangerspctrm(row,col,:)),'color',cmap(ncov,:),'LineWidth',2)
            hold on
            if col == 1
                ylabel(['signal ' num2str(row)])
            end
            if row ==Nsig
                xlabel(['signal ' num2str(col)])
            end
            if row == 2 && col ==1
                ylabel({'from:' ; ['signal ' num2str(row)]})
            end
            ylim([0 1])
        end
    end
    if ncov == NC
        set(gcf,'Position',[636 172 1264 889])
        h = legend(NCtits);
        set(h,'Position',[0.895 0.845 0.0837 0.0804],'FontSize',8)
        annotation(gcf,'textbox',...
            [0.351 0.941 0.321 0.043],...
            'String',{'NonParametric Granger'},...
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
            x = []; y = []; z = []; z2 = []
            for p = 1:size(data.trial,2)
                x = [x data.trial{p}(i,:)];
                y = [y data.trial{p}(j,:)];
                z = [z data.trial{p}(2,:)];
                z2= [z2 data.trial{p}(2:3,:)];
            end
            %                 [f13,~,~]=sp2a2_R2_mt(x',y',FTdata.fsample,7,'M1');
            [f13,t13,c13]=sp2_R2a_pc1(x,y,z,data.fsample,2^8);
            npdspctrmZ{1,1}(i,j,:) = f13(:,10);
            npdspctrmZ{1,2}(i,j,:) = f13(:,12); % Backward
            npdspctrmZ{1,3}(i,j,:) = f13(:,11); % Forward
            
            [f1z2,t1z2] = HOpartcohtw_160517(x',y',z2',250,8,0);
            npdspctrmZ2(i,j,:)= f1z2(:,4);
            
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
            if i ==Nsig
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
        set(gcf,'Position',[636 172 1264 889])
        h = legend(NCtits);
        set(h,'Position',[0.895 0.845 0.0837 0.0804],'FontSize',8)
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
            if i ==Nsig
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
        set(gcf,'Position',[636 172 1264 889])
        h = legend(NCtits);
        set(h,'Position',[0.895 0.845 0.0837 0.0804],'FontSize',8)
        annotation(gcf,'textbox',...
            [0.351 0.941 0.321 0.043],...
            'String',{'NonParametric Directionality XCOR'},...
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
            if i ==Nsig
                xlabel(['signal ' num2str(j)])
            end
            if i == 2 && j ==1
                ylabel({'from:' ; ['signal ' num2str(i)]})
            end
            if ncov == NC
                plot(fx,repmat(c13.ch_c95,1,size(fy,1)),'k--')
            end
            ylim([0 1])
            xlim([0 100])
            clear fy
        end
    end
    if ncov==NC
        set(gcf,'Position',[636 172 1264 889])
        h = legend(NCtits);
        set(h,'Position',[0.895 0.845 0.0837 0.0804],'FontSize',8)
        annotation(gcf,'textbox',...
            [0.351 0.941 0.321 0.043],...
            'String',{'NonParametric Directionality Partialed on 1: Spectra'},...
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
            if i ==Nsig
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
        set(gcf,'Position',[636 172 1264 889])
        h = legend(NCtits);
        set(h,'Position',[0.895 0.845 0.0837 0.0804],'FontSize',8)
        annotation(gcf,'textbox',...
            [0.351 0.941 0.321 0.043],...
            'String',{'NonParametric Directionality RhoYX Partialed on 1: XCOR'},...
            'HorizontalAlignment','center',...
            'FitBoxToText','off',...
            'EdgeColor','none');
    end
    %% NPDZ2 Spectrum
    figure(9)
    for i = 1:length(data.label)
        for j = 1:length(data.label)
            linind = sub2ind([length(data.label),length(data.label)],j,i);
            subplot(length(data.label),length(data.label),linind)
            fy = squeeze(npdspctrmZ2(i,j,:));
            fx = f1z2(:,1);
            
            plot(fx,fy,'color',cmap(ncov,:),'LineWidth',2)
            hold on
            if j == 1
                ylabel(['signal ' num2str(i)])
            end
            if i ==Nsig
                xlabel(['signal ' num2str(j)])
            end
            if i == 2 && j ==1
                ylabel({'from:' ; ['signal ' num2str(i)]})
            end
            %             if ncov == NC
            %                 plot(fx,repmat(c13.rho_c95,1,size(fy,1)),'k--')
            %                 hold on
            %                 plot(fx,repmat(-c13.rho_c95,1,size(fy,1)),'k--')
            %             end
            ylim([0 1])
            xlim([0 100])
            clear fy
        end
    end
    if ncov==NC
        set(gcf,'Position',[636 172 1264 889])
        h = legend(NCtits);
        set(h,'Position',[0.895 0.845 0.0837 0.0804],'FontSize',8)
        annotation(gcf,'textbox',...
            [0.351 0.941 0.321 0.043],...
            'String',{'2nd Order Partialised Coherence conditioned on 1 and 2'},...
            'HorizontalAlignment','center',...
            'FitBoxToText','off',...
            'EdgeColor','none');
    end
end
saveallfiguresFIL('C:\Users\twest\Documents\Work\PhD\NPD_Validate\MVAR_model\figures\4link\outputs_4link','-jpg')