function [Hz granger grangerft] = computeGranger(freq,cmap,Nsig,plotfig,linestyle,mvflag,bstrap)
if nargin<6
    mvflag = 0;
end
if nargin<7
    bstrap = 0;
end
null = nan(size(freq.freq));
if mvflag == 1
    %Multivariate
    [grangerft,igrangerft] = NPG_freq(freq,[],bstrap);
end
for i = 1:length(freq.label)
    for j = 1:length(freq.label)
        if mvflag == 0 % Pairwise
            if i == j
                %                 [grangerft,igrangerft] = NPG_freq(freq,[i j],0);
            else
                [grangerft,igrangerft] = NPG_freq(freq,[i j],bstrap);
            end
            if i<j
                ip = 1; jp = 2;
            else
                ip = 2; jp = 1;
            end
        else
            ip = i; jp = j;
        end
        if i==j
            granger{1,1}(i,j,:) = null; %squeeze(igrangerft.instantspctrm(1,1,:));
            granger{1,2}(i,j,:) = null; %squeeze(grangerft.grangerspctrm(1,1,:)); % Backward (j -> i)
            granger{1,3}(i,j,:) = null; %squeeze(grangerft.grangerspctrm(1,1,:)); % Forward (i -> j)
            
            granger{2,1}(i,j,:) = null; %squeeze(igrangerft.cispec(1,1,:));
            granger{2,2}(i,j,:) = null; %squeeze(grangerft.cispec(1,1,:)); % Backward (j -> i)
            granger{2,3}(i,j,:) = null; %squeeze(grangerft.cispec(1,1,:)); % Forward (i -> j)
        else
            granger{1,1}(i,j,:) = squeeze(igrangerft.instantspctrm(ip,jp,:));
            granger{1,2}(i,j,:) = squeeze(grangerft.grangerspctrm(jp,ip,:)); % Backward (j -> i)
            granger{1,3}(i,j,:) = squeeze(grangerft.grangerspctrm(ip,jp,:)); % Forward (i -> j)
            
            granger{2,1}(i,j,:) = squeeze(igrangerft.cispec(ip,jp,:));
            granger{2,2}(i,j,:) = squeeze(grangerft.cispec(jp,ip,:)); % Backward (j -> i)
            granger{2,3}(i,j,:) = squeeze(grangerft.cispec(ip,jp,:)); % Forward (i -> j)
        end
    end
end


if plotfig
%     plotNPD(grangerft.freq,granger,freq,cmap,1,linestyle,bstrap)
    plotDiffNPD(grangerft.freq,granger,freq,cmap,1,linestyle,bstrap)
end


% Put in same format as the NPD
% granger{1,1} = igrangerft.instantspctrm;
% granger{1,2} = grangerft.grangerspctrm(2:-1:1,2:-1:1,:);
% granger{1,3} = grangerft.grangerspctrm;
Hz = grangerft.freq;
% % % Plotting if multivariate
% %     for row= 1:Nsig
% %         for col=1:Nsig
% %             subplot(Nsig,Nsig,(row-1)*Nsig+col);
% %             plot(grangerft.freq, squeeze(grangerft.grangerspctrm(row,col,:)),'color',cmap,'LineWidth',2,'linestyle',linestyle)
% %             hold on
% %             if col == 1
% %                 ylabel(['signal ' num2str(row)])
% %             end
% %             if row == Nsig
% %                 xlabel(['signal ' num2str(col)])
% %             end
% %             if row == Nsig && col == floor(Nsig/2)
% %                 xlabel({'to:' ; ['signal ' num2str(col)]})
% %             end
% %             if row == floor(Nsig/2) && col == 1
% %                 ylabel({'from:' ; ['signal ' num2str(row)]})
% %             end
% %             ylim([0 1])
% %         end
% %     end
