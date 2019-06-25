function score = matrixScore(Bc,Z)
% B is the predicted; Z is the actual; crit is the
SB = Bc+Z;
score = 100.*(sum(SB(:)== 2 | SB(:)== 0) - sum(SB(:)== 1) - size(diag(Bc),1) )./(numel(Bc) - size(diag(Bc),1));