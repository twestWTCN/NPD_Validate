close all
% Results of confidence 
% Type I - Shuffling FFTs
A = [   -1.0000    0.0066    0.0047
   -1.3010    0.0091    0.0062
   -2.0000    0.0151    0.0095
   -3.0000    0.0240    0.0145];

plot(A(:,1),A(:,2),'r-O'); % NPD
hold on
plot(A(:,1),A(:,3),'b-O'); % NPG

% Type II - Phase Randomization
A = [   -1.0000    0.0178    0.0120
   -1.3010    0.0228    0.0148
   -2.0000    0.0345    0.0214
   -3.0000    0.0512    0.0304];

plot(A(:,1),A(:,2),'r--S');
hold on
plot(A(:,1),A(:,3),'b--S')
grid on
xlabel('log_{10} Significance Level'); ylabel('Rank Order Significance Threshold')
legend({'FFT Reshuffle - NPD','FFT Reshuffle - NPG','Phase Randomization - NPD','Phase Randomization - NPG'})