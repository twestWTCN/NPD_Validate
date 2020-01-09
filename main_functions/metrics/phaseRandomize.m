function dxShuff = phaseRandomize(dx)
% takes in a L samples x N trials Fourier transform and randomizes the
% phase whilst preserving magnitude
[L,N] = size(dx);
if rem(L,2)
    L = L - 1;
end
randPhase = randbetween(0,2*pi,L/2-1,N);
dx(2:L/2,:) = dx(2:L/2,:).*exp(1i*randPhase);
dx(L/2+2:L,:) = dx(L/2+2:L,:).*exp(-1i*flipud(randPhase));
dxShuff = dx;