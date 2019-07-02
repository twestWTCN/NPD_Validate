function [f13 t13 f13Z t13Z f13W t13W f13Q t13Q ci13 ci13Z ci13W ci13Q] = NPD_XYZWQ(x,y,z,w,q,fsamp,npdord,bstrp)
%                 [f13,~,~]=sp2a2_R2_mt(x',y',FTdata.fsample,7,'M1');
[f13,t13,~]=sp2a2_R2_tw(x',y',fsamp,npdord);
winsize = fix(2^npdord);
if rem(winsize,2)
    winsize = winsize +1;
end
[f13Z,t13Z,~]=sp2_R2a_pc1_tw(x',y',z',fsamp,winsize);
[f13W,t13W,~]=sp2_R2a_pc1_tw(x',y',w',fsamp,winsize);
[f13Q,t13Q,~]=sp2_R2a_pc1_tw(x',y',q',fsamp,winsize);

if bstrp == 1
    bsn = 100; %200;
    disp('bsn is v low!')
    sf13 = nan([bsn,size(f13)]);
    st13 = nan([bsn,size(t13)]);
    sf13Z= nan([bsn,size(f13Z)]);
    st13Z= nan([bsn,size(t13Z)]);
    sf13W= nan([bsn,size(f13W)]);
    st13W= nan([bsn,size(t13W)]);
    sf13Q= nan([bsn,size(f13Q)]);
    st13Q= nan([bsn,size(t13Q)]);    
    parfor i = 1:bsn
        xshuff = vectorShuff(x,winsize);
        yshuff = vectorShuff(y,winsize);
        zshuff = vectorShuff(z,winsize);
        wshuff = vectorShuff(w,winsize);
        [sf13(i,:,:),st13(i,:,:),sf13Z(i,:,:),st13Z(i,:,:),sf13W(i,:,:),st13W(i,:,:),st13Q(i,:,:),sf13Q(i,:,:)] = NPD_XYZWQ(xshuff,yshuff,zshuff,wshuff,qshuff,fsamp,npdord,0);
    end
    ci13 = repmat(mean(squeeze(prctile(sf13,99.9,1)),1),size(sf13,2),1);
    ci13Z = repmat(mean(squeeze(prctile(sf13Z,99.9,1)),1),size(sf13,2),1);
    ci13W = repmat(mean(squeeze(prctile(sf13W,99.9,1)),1),size(sf13,2),1);
    ci13Q = repmat(mean(squeeze(prctile(sf13Q,99.9,1)),1),size(sf13,2),1);
else
    ci13 = zeros(size(f13));
    ci13Z = zeros(size(f13Z));
    ci13W = zeros(size(f13W));
    ci13Q = zeros(size(f13Q));
end