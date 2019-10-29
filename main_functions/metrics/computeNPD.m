function [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv npdcrcvZ npdcrcvW] = computeNPDQ(data,frstord,npdord,sclr,bstrap)
if nargin<3
    npdord = 8;
end
if nargin<4
    sclr = 1;
end
if nargin<5
    bstrap = 0;
end
% Partialised
for i = 1:length(data.label)
    for j = 1:length(data.label)
        x = []; y = []; z = []; w = [];
        for p = 1:size(data.trial,2)
            x = [x data.trial{p}(i,:)];
            y = [y data.trial{p}(j,:)];
            z = [z data.trial{p}(frstord,:)];
            w = [w data.trial{p}(frstord+1,:)];
%             q = [w data.trial{p}(frstord+2,:)];
            %                 z2= [z2 data.trial{p}(sndord,:)];
        end
        fsamp = data.fsample;
        [f13 t13 f13Z t13Z f13W t13W ci13 ci13Z ci13W] = NPD_XYZW(x,y,z,w,fsamp,npdord,bstrap);
%         [f13 t13 f13Z t13Z f13W t13W f13Q t13Q ci13 ci13Z ci13W ci13Q] = NPD_XYZWQ(x,y,z,w,q,fsamp,npdord,bstrap);
        npdspctrm{1,1}(i,j,:) = sclr*f13(:,10);
        npdspctrm{1,2}(i,j,:) = sclr*f13(:,12); % Backward (j -> i)
        npdspctrm{1,3}(i,j,:) = sclr*f13(:,11); % Forward (i -> j)
        
        npdspctrm{2,1}(i,j,:) = sclr*ci13(:,10);
        npdspctrm{2,2}(i,j,:) = sclr*ci13(:,12);
        npdspctrm{2,3}(i,j,:) = sclr*ci13(:,11);
        
        npdcrcv(i,j,:) = t13(:,3);
        nscohspctrm{1}(i,j,:) = f13(:,4);
        nscohspctrm{2}(i,j,:) = ci13(:,4);
        
        npdspctrmZ{1,1}(i,j,:) = sclr*f13Z(:,10);
        npdspctrmZ{1,2}(i,j,:) = sclr*f13Z(:,12); % Backward
        npdspctrmZ{1,3}(i,j,:) = sclr*f13(:,11); % Forward
        npdcrcvZ(i,j,:) = t13(:,3);
        
        npdspctrmZ{2,1}(i,j,:) = sclr*ci13Z(:,10);
        npdspctrmZ{2,2}(i,j,:) = sclr*ci13Z(:,12);
        npdspctrmZ{2,3}(i,j,:) = sclr*ci13Z(:,11);        
        
        npdspctrmW{1,1}(i,j,:) = sclr*f13W(:,10);
        npdspctrmW{1,2}(i,j,:) = sclr*f13W(:,12); % Backward
        npdspctrmW{1,3}(i,j,:) = sclr*f13W(:,11); % Forward
        npdcrcvW(i,j,:) = t13W(:,3);
        
        npdspctrmW{2,1}(i,j,:) = sclr*ci13W(:,10);
        npdspctrmW{2,2}(i,j,:) = sclr*ci13W(:,12);
        npdspctrmW{2,3}(i,j,:) = sclr*ci13W(:,11);        
   
%         npdspctrmQ{1,1}(i,j,:) = sclr*f13Q(:,10);
%         npdspctrmQ{1,2}(i,j,:) = sclr*f13Q(:,12); % Backward
%         npdspctrmQ{1,3}(i,j,:) = sclr*f13Q(:,11); % Forward
%         npdcrcvQ(i,j,:) = t13Q(:,3);
%            
%         
        
    end
end
Hz = f13(:,1);
lags = t13(:,1);
