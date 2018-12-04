function [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv npdcrcvZ npdcrcvW] = computeNPD(data,frstord,npdord)
if nargin<3
    npdord = 8;
end
% Partialised
for i = 1:length(data.label)
    for j = 1:length(data.label)
        %                 ftspect = freqFour.fourierspctrm(:,[chcombs(1,i) chcombs(2,i)],:);
        %                 [f13,t13,cl13] = sp2a2_R2_TW(FTdata.fsample,8,ftspect);
        x = []; y = []; z = []; w = [];
        for p = 1:size(data.trial,2)
            x = [x data.trial{p}(i,:)];
            y = [y data.trial{p}(j,:)];
            z = [z data.trial{p}(frstord,:)];
            w = [w data.trial{p}(frstord+1,:)];
            %                 z2= [z2 data.trial{p}(sndord,:)];
        end
        %                 [f13,~,~]=sp2a2_R2_mt(x',y',FTdata.fsample,7,'M1');
        [f13,t13,c13]=sp2a2_R2(x',y',data.fsample,npdord);
        npdspctrm{1,1}(i,j,:) = f13(:,10);
        npdspctrm{1,2}(i,j,:) = f13(:,12); % Backward
        npdspctrm{1,3}(i, j,:) = f13(:,11); % Forward
        npdcrcv(i,j,:) = t13(:,3);
        nscohspctrm(i,j,:) = f13(:,4);
        
        [f13,t13,c13]=sp2_R2a_pc1(x',y',z',data.fsample,2^npdord);
        npdspctrmZ{1,1}(i,j,:) = f13(:,10);
        npdspctrmZ{1,2}(i,j,:) = f13(:,12); % Backward
        npdspctrmZ{1,3}(i,j,:) = f13(:,11); % Forward
        npdcrcvZ(i,j,:) = t13(:,3);
        
        [f13,t13,c13]=sp2_R2a_pc1(x',y',w',data.fsample,2^npdord);
        npdspctrmW{1,1}(i,j,:) = f13(:,10);
        npdspctrmW{1,2}(i,j,:) = f13(:,12); % Backward
        npdspctrmW{1,3}(i,j,:) = f13(:,11); % Forward
        npdcrcvW(i,j,:) = t13(:,3);
        %             [f1z2,t1z2] = HOpartcohtw_160517(x',y',z2',250,8,0);
        %             npdspctrmZ2(i,j,:)= f1z2(:,4);
        
        %             npdcrcvZ(i,j,:) = t13(:,3);
        %             nscohspctrmZ(i,j,:) = f13(:,4);
    end
end
Hz = f13(:,1);
lags = t13(:,1);
