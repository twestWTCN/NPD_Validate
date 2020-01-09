function [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv npdcrcvZ npdcrcvW] = ft_computeNPD(freq,fsamp,frstord,npdord,perm,permtype)
%% Function to compute NPD using modified Neurospec Scripts taking Fieldtrip Freq Structure
if nargin<3
    npdord = 8;
end
if nargin<5
    perm = 0;
end
if nargin<6
    permtype = 0;
end
% Partialised
for i = 1:length(freq.label)
    for j = 1:length(freq.label)
        dx = squeeze(freq.fourierspctrm(:,i,:));
        dy = squeeze(freq.fourierspctrm(:,j,:));
        dz = squeeze(freq.fourierspctrm(:,frstord,:));
        dw = squeeze(freq.fourierspctrm(:,frstord+1,:));
        
        [f13 t13 f13Z t13Z f13W t13W ci13 ci13Z ci13W] = ft_NPD_XYZW(dx,dy,dz,dw,fsamp,npdord,perm,permtype);
        
        npdspctrm{1,1}(i,j,:) = f13(:,10);
        npdspctrm{1,2}(i,j,:) = f13(:,11); % Forward (i -> j)
        npdspctrm{1,3}(i,j,:) = f13(:,12); % Backward (j -> i)
        
        npdspctrm{2,1}(i,j,:) = ci13(:,10);
        npdspctrm{2,2}(i,j,:) = ci13(:,11);
        npdspctrm{2,3}(i,j,:) = ci13(:,12);
        
        npdcrcv(i,j,:) = t13(:,3);
        nscohspctrm{1}(i,j,:) = f13(:,4);
        nscohspctrm{2}(i,j,:) = ci13(:,4);
        
        npdspctrmZ{1,1}(i,j,:) = f13Z(:,10);
        npdspctrmZ{1,2}(i,j,:) = f13Z(:,11); 
        npdspctrmZ{1,3}(i,j,:) = f13(:,12); 
        npdcrcvZ(i,j,:) = t13(:,3);
        
        npdspctrmZ{2,1}(i,j,:) = ci13Z(:,10);
        npdspctrmZ{2,2}(i,j,:) = ci13Z(:,11);
        npdspctrmZ{2,3}(i,j,:) = ci13Z(:,12);
        
        npdspctrmW{1,1}(i,j,:) = f13W(:,10);
        npdspctrmW{1,2}(i,j,:) = f13W(:,11); % Backward
        npdspctrmW{1,3}(i,j,:) = f13W(:,12); % Forward
        npdcrcvW(i,j,:) = t13W(:,3);
        
        npdspctrmW{2,1}(i,j,:) = ci13W(:,10);
        npdspctrmW{2,2}(i,j,:) = ci13W(:,11);
        npdspctrmW{2,3}(i,j,:) = ci13W(:,12);
        
        %         npdspctrmQ{1,1}(i,j,:) = sclr*f13Q(:,10);
        %         npdspctrmQ{1,2}(i,j,:) = sclr*f13Q(:,12); % Backward
        %         npdspctrmQ{1,3}(i,j,:) = sclr*f13Q(:,11); % Forward
        %         npdcrcvQ(i,j,:) = t13Q(:,3);
    end
end
Hz = f13(:,1);
lags = t13(:,1);
