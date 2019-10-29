function [b,a] = buildfilter()
[n,Wn]  = buttord([30 60]/200, [25 65]/200, 1, 10);
[b,a]   = butter(n,Wn,'bandpass');
