function [noise, targ, filter, SNR, delay, filters] = find_in_mix(dmix, dclean, srmix, Tfilt, Tpre, Twin, Thop)
% [noise, targ, filter, SNR, delay] = find_in_mix(dmix, dclean, srmix, Tfilt, Tpre, Twin, Thop)
%     dmix consists of a filtered version of dclean with added
%     noise.  Use best linear fit to return the noise without the
%     clean signal, the filtered (and slightly time-warped) version of clean
%     that best fits it, and the filter that can be applied to the clean
%     signal to make it match what was in the mix, and the overall
%     SNR this implies for the mixture.
%     Tfilt is the duration of the FIR filter to estimate (0.040 s), 
%     and Tpre is how much "pre-echo" to allow before the timing
%     peak (0.005 s).
%     Twin is the length of windows used for estimating filter
%     (default 8 s), and Thop is the hop between successive windows
%     (default Twin/2).
% 2011-02-10 Dan Ellis dpwe@ee.columbia.edu

if nargin < 4
  Tfilt = 0.040; % duration of FIR filter
  %Tfilt = 0.100; % duration of FIR filter
end

Lfilt = round(Tfilt*srmix);
if nargin < 5
  Tpre  = 0.005; % how much to allow before peak
  %Tpre  = 0.020; % how much to allow before peak
end

if nargin < 6; Twin = 8.0; end
if nargin < 7; Thop = Twin/2; end

% find best skew (up to half a second)
skewmaxsec = 0.5;
skewsamp = round(skewmaxsec*srmix);
%resolution = 1;
% so as not to hit limit in fftfilt.. under MATLAB R2010b
resolution = ceil(length(dclean)/2^20);
dosquare = 0;
mixdelay = find_skew(dmix, dclean, skewsamp, resolution, dosquare);

delay = mixdelay/srmix;
disp(['Delay = ',sprintf('%.6f',delay),' s']);

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

%Twin = 1.0;  % match on 1 sec blocks
%Twin = 8.0;
% longer gives better results for signals with no time warp issues
Lwin = round(Twin*srmix);
%Thop = Twin/2;
Lhop = round(Thop*srmix);

[targ, noise, filters, Es] = decomp_lin_win(dmix, dclean, Lfilt, ...
                                            Lwin, Lhop);

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