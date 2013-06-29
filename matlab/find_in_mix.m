function [noise, targ, filter, SNR, delay, filters] = find_in_mix(dmix, ...
                                                  dclean, srmix, Tfilt)
% [noise, targ, filter, SNR, delay] = find_in_mix(dmix, dclean, srmix)
%     dmix consists of a filtered version of dclean with added
%     noise.  Use best linear fit to return the noise without the
%     clean signal, the filtered (and slightly time-warped) version of clean
%     that best fits it, and the filter that can be applied to the clean
%     signal to make it match what was in the mix, and the overall
%     SNR this implies for the mixture.
% 2011-02-10 Dan Ellis dpwe@ee.columbia.edu

% find best skew (up to half a second)
skewmaxsec = 15.0;
%skewmaxsec = 0.5;
mixdelay = find_skew(dclean, dmix, round(skewmaxsec*srmix));

delay = mixdelay/srmix;
disp(['Delay = ',num2str(delay),' s']);

if nargin < 4
  Tfilt = 0.040; % duration of FIR filter
  %Tfilt = 0.100; % duration of FIR filter
end

Lfilt = round(Tfilt*srmix);
Tpre  = 0.005; % how much to allow before peak
%Tpre  = 0.020; % how much to allow before peak
Lpre  = round(Tpre*srmix);

%disp(['skewmaxsec=',num2str(skewmaxsec),' Tfilt=',num2str(Tfilt)]);

skew = mixdelay - Lpre;

% align signals
if skew > 0
  % mix is delayed relative to clean, so chop off its start
  dmix = dmix((skew+1):end);
else
  % clean is actually delayed relative to mix, so chop off its start
  dclean = dclean((1-skew):end);
end

% % make files same length
% dlen = min(length(dmix),length(dclean));
% dmax = dmax(1:dlen);
% dclean = dclean(1:dlen);

%Twin = 1.0;  % match on 2 sec blocks
Twin = 0.5;  % match on 2 sec blocks
Lwin = round(Twin*srmix);
[targ, noise, filters, Es] = decomp_lin_win(dmix, dclean, Lfilt, Lwin);

%% Choose just one filter to return; the "median"
%fnorms = sum(filters.^2);
%% sort the norms and take the middle one
%[vv,xx] = sort(fnorms);
%bestf = xx(round(length(xx)/2));

% No, take the filter from the frame where the excitation had the
% most energy
[biggestE, bestf] = max(Es);

filter = filters(:,bestf);

% figure the SNR by finding the "active level" of signal and noise
siglevel = activlev(targ, srmix);
noiselevel = activlev(noise, srmix);
SNR = 10 * log10(siglevel/noiselevel);

disp(['SNR = ', num2str(SNR),' dB']);