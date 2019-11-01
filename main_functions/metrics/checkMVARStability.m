function [eigVals,flag] = checkMVARStability(C)
flag = 0;
mat = []; row = [];
for pout = 1:size(C,3)
    for pin = 1:size(C,3)
        if pout == 1
            row = [row squeeze(C(:,:,pin))];
        elseif (pout-1) == pin
            row = [row eye(size(C,1))];
        else
            row = [row zeros(size(C,1))];
        end
    end
    mat = [mat;row]; row = [];
end
CM = mat;
eigVals = abs(eig(CM));
if any(eigVals>1)
    flag = 1;
    disp('Model is not stable!')
end

% CM = [squeeze(C(:,:,1)) squeeze(C(:,:,2)) squeeze(C(:,:,3));
%     
% eye(size(C,1))    zeros(size(C,1))    zeros(size(C,1));
% 
% zeros(size(C,1))  eye(size(C,1))       zeros(size(C,1))];
% 



