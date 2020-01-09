function [f13 t13 f13Z t13Z f13W t13W ci13 ci13Z ci13W] = ft_NPD_XYZW(dx,dy,dz,dw,fsamp,npdord,bstrp,bstrptype)
%                 [f13,~,~]=sp2a2_R2_mt(x',y',FTdata.fsample,7,'M1');
[f13,t13,~]=ft_sp2a2_R2_tw(dx',dy',fsamp,npdord,0);

[f13Z,t13Z,~]=ft_sp2_R2a_pc1_tw(dx',dy',dz',fsamp,npdord,0);
[f13W,t13W,~]=ft_sp2_R2a_pc1_tw(dx',dy',dw',fsamp,npdord,0);

if bstrp == 1
    bsn = 100; %200;
    disp('bsn is v low!')
    sf13 = nan([bsn,size(f13)]);
    st13 = nan([bsn,size(t13)]);
    sf13Z= nan([bsn,size(f13Z)]);
    st13Z= nan([bsn,size(t13Z)]);
    sf13W= nan([bsn,size(f13W)]);
    st13W= nan([bsn,size(t13W)]);
    parfor i = 1:bsn
        [sf13(i,:,:),st13(i,:,:),~]=ft_sp2a2_R2_tw(dx',dy',fsamp,npdord,bstrptype);
        [sf13Z(i,:,:),st13Z(i,:,:),~]=ft_sp2_R2a_pc1_tw(dx',dy',dz',fsamp,npdord,bstrptype);
        [sf13W(i,:,:),st13W(i,:,:),~]=ft_sp2_R2a_pc1_tw(dx',dy',dw',fsamp,npdord,bstrptype);
    end
    ci13 = repmat(mean(squeeze(prctile(sf13,99.9,1)),1),size(sf13,2),1);
    ci13Z = repmat(mean(squeeze(prctile(sf13Z,99.9,1)),1),size(sf13,2),1);
    ci13W = repmat(mean(squeeze(prctile(sf13W,99.9,1)),1),size(sf13,2),1);
else
    ci13 = zeros(size(f13));
    ci13Z = zeros(size(f13Z));
    ci13W = zeros(size(f13W));
end