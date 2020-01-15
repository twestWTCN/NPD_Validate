function [] = wrapper_Fig5A_SymSNR(C,NCV,NC,permrun,fresh)
% Sweep across symmetrical SNRs
NCvec = logspace(log10(0.001),log10(1000),NC);
NCvec = sqrt(NCvec); % Convert to std
[~,NCbase] = min(abs(1-NCvec)); % Find base condition
cmap = linspecer(4);

% Run the MVAR Simulation (fieldtrip implementation uses BSMART)
rng(12312)
Nsig = size(C,1);
cfg             = [];
cfg.ntrials     = 1;
cfg.triallength = 250;
cfg.fsample     = 200;
cfg.nsignal     = Nsig;
cfg.method      = 'ar';
cfg.params = C;
cfg.noisecov = NCV;
X              = ft_connectivitysimulation(cfg);
segOrd = 8; % 2^n length of segment used for FFT

for ncov = 1:NC
        disp(ncov)

    if ncov == NCbase
        if fresh == 1
            perm = 1;
            permtype = permrun;
        elseif fresh == 0
            % Can switch off if you have run once to save the results see
            % L75-76
            perm = 0;
            permtype = 0;
        end
    else
        perm = 0;
        permtype = 0;
    end
      
    % Plotting parameters
    linestyle = '-';
    cmapn = cmap;
    plotfig =0;
    
    %% Compute Symmetrical SNR
    data = X;
    ncv = NCvec(ncov);
    
    randproc = randn(size(data.trial{1}));
    for i = 1:size(randproc,1)
        s = data.trial{1}(i,:);
        s = (s-mean(s))./std(s);
        n = (ncv.*randproc(i,:));
        y = s + n;
        snr = var(s)/var(n);
        snrbank(ncov,i) = snr;
        data.trial{1}(i,:) = y;
    end
    
    %% Compute Power with Fieldtrip (used for NPD and NPG computations)
    datalength = (2^segOrd)./cfg.fsample;
    freq = computeSpectra(data,[0 0 0],Nsig,plotfig,linestyle,1,datalength);
    pow = mean(abs(squeeze(freq.fourierspctrm(:,1,:))),1); % Convert to power
    npPow(ncov) = max(pow(freq.freq>42 & freq.freq<62));

    %% Compute NPD with Neurospec
    coh = [];
    [Hz lags npdspctrm npdspctrmZ npdspctrmW nscohspctrm npdcrcv] = ft_computeNPD(freq,cfg.fsample,1,segOrd,perm,permtype);
    coh.freq= Hz; coh.cohspctrm = nscohspctrm{1};
    nscoh(ncov) = max(nscohspctrm{1}(1,2,Hz>42 & Hz<62));
    npd(ncov) = max(npdspctrm{1,3}(1,2,Hz>42 & Hz<62));
    npdCi(ncov) = mean(npdspctrm{2,2}(1,2,:));

    %% Compute NPG
    [Hz granger grangerft] = computeGranger(freq,0,perm,permtype);
    npGC(ncov) = max(granger{1,3}(1,2,grangerft.freq>42 & grangerft.freq<62));
    npGCci(ncov) = mean(granger{2,3}(1,2,:));
end

if fresh == 1
    % save(['C:\Users\Tim\Documents\Work\GIT\NPD_Validate\precomp_CI_table\F5A_CItab_type' num2str(permrun)],'npdCi','npGCci')
    save(['C:\Users\timot\Documents\GitHub\NPD_Validate\precomp_CI_table\F5A_CItab_type' num2str(permrun)],'npdCi','npGCci')
elseif fresh == 0
    % load(['C:\Users\Tim\Documents\Work\GIT\NPD_Validate\precomp_CI_table\F5A_CItab_type' num2str(permrun)],'npdCi','npGCci')
    load(['C:\Users\timot\Documents\GitHub\NPD_Validate\precomp_CI_table\F5A_CItab_type' num2str(permrun)],'npdCi','npGCci')
end
SNRDB = 10*log10(snrbank(:,1));

A = SNRDB;
% scatter(NCvec,npPow,40,cmapn(1,:),'filled')
pa(1) = scatter(A,nscoh,50,cmapn(1,:),'Marker','+','LineWidth',3);;
hold on
% [param,stat]=sigm_fit(A,nscoh)
% hold on
pa(2) = scatter(A,npd,40,cmapn(3,:),'filled');
plot(A,repmat(npdCi(NCbase),1,size(A,1)),'LineStyle','--','color',cmapn(3,:))

% [param,stat]=sigm_fit(A,npd)
pa(3) = scatter(A,npGC,40,cmapn(2,:),'filled');
plot(A,repmat(npGCci(NCbase),1,size(A,1)),'LineStyle','--','color',cmapn(2,:))

% [param,stat]=sigm_fit(A,npGC); %,[],[.0043 0.5884 0.689 -3.76])
grid on
xlabel('SNR_{dB}');ylabel('FC Magnitude')
legend(pa,{'Coherence','NPD','Granger'})
title('FC vs SNR')
xlim([-inf inf])