function [grangerft,igrangerft] = NPG_freq(freq,chinds,bstrap,bstraptype)
cfg           = [];
if ~isempty(chinds)
    cfg.channel   = freq.label(chinds); 
end
cfg.method    = 'granger';
grangerft       = ft_connectivityanalysis(cfg, freq);

cfg           = [];
if ~isempty(chinds)
    cfg.channel   = freq.label(chinds);
end
cfg.method    = 'instantaneous_causality';
igrangerft      = ft_connectivityanalysis(cfg, freq);

if bstrap == 1 && bstraptype>0
    bsn = 1000;
    ispec = nan([size(igrangerft.instantspctrm) bsn]);
    spec = nan([size(grangerft.grangerspctrm) bsn]);
    
    parfor i = 1:bsn
        freqS = freqShuff(freq,bstraptype); % If bootstrapping spectra
%         freqS = phaseRand(freq,bstraptype); % If phase randomizing
        [grangerft,igrangerft] = NPG_freq(freqS,chinds,0);
        ispec(:,:,:,i) = igrangerft.instantspctrm;
        spec(:,:,:,i) = grangerft.grangerspctrm;
    end
    igrangerft.cispec = meanVec(ispec);
    grangerft.cispec =  meanVec(spec);
else
    igrangerft.cispec = zeros(size(igrangerft.instantspctrm));
    grangerft.cispec = zeros(size(grangerft.grangerspctrm));
end

function L = meanVec(spec)
for i = 1:size(spec,1);
    for j = 1:size(spec,2)
        L(i,j,:) = repmat(mean(squeeze(prctile(spec(i,j,:,:),99.99,4))),1,size(spec,3));
    end
end