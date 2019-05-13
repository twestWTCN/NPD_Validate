function C = placeCon(C,csel,osel,cstrength)
for i = 1:size(csel,2)
    C(csel(1,i),csel(2,i),osel(i)) = cstrength;
end