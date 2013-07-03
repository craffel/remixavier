function [D, SR] = audioread(FN, TARGETSR, FORCEMONO, START, DUR)
% [D, SR] = audioread(FN, TARGETSR, FORCEMONO, START, DUR)
%   Read in an audio file, using wavread, mp3read, m4aread, or flacread as
%   appropriate. 
% 2010-09-16 Dan Ellis dpwe@ee.columbia.edu
% $Header: /Users/dpwe/matlab/columbiafns/RCS/audioread.m,v 1.2 2011/09/09 04:13:18 dpwe Exp dpwe $

if nargin < 2; TARGETSR = []; end
if nargin < 3; FORCEMONO = 0; end
if nargin < 4; START = 0; end
if nargin < 5; DUR = 0; end

END = START+DUR;

if length(TARGETSR) == 0
  TARGETSR = 0;
end

if length(FN) > 2 && FN(1) == '~' && FN(2) == '/'
  % tilde substitution
  FN = [getenv('HOME'),FN(2:end)];
end


if exist(FN, 'file') == 0
  error(['audioread: file ',FN,' not found']);
end

[pth,nam,ext] = fileparts(FN);
ext = lower(ext);

ismp3 = strcmp(ext,'.mp3');
ism4a = strcmp(ext, '.m4a') || strcmp(ext, '.aac') || strcmp(ext, '.mp4');
iswav = strcmp(ext, '.wav');
isrds = strcmp(ext, '.wv1') || strcmp(ext, '.wv2') ...
        || strcmp(ext, '.wvs') || strcmp(ext, '.shn') ...
        || strcmp(ext, '.aif') || strcmp(ext, '.aiff') ...
        || strcmp(ext, '.sph');
isflac = strcmp(ext, '.flac');

if ismp3 || ism4a || iswav || isrds || isflac
  if ismp3
    [SIZ,SR] = mp3read(FN,'size');
  elseif ism4a
    [SIZ,SR] = m4aread(FN,'size');
  elseif iswav
    [SIZ,SR] = wavread(FN,'size');
  elseif isflac
    [SIZ,SR] = flacread(FN,'size');
  else
    % ReadSound has no 'size' function
    SIZ = [-1 0];
    % but we really want to know the SR - read one sample
    [D,SR] = ReadSound(FN,1);
  end

  ds = 1;
  % apply downsampling during read to minimize core usage
  % (not supported in ReadSound)
  if ~isrds && TARGETSR > 0 && SR >= 2*TARGETSR
    ds = 2;
    if SR >= 4*TARGETSR && (SR >= 22050 || ~ismp3)
      % mpg123 will only do 4:1 downsampling for higher SRs
      ds = 4;
    end
  end
  % apply truncation during read to avoid massive core usage
%  N = 0;
%  if START > 0 || END > START
    SIX = max(1,round(START*SR/ds)+1);
    if ~ismp3 && ~ism4a && END <= START && SIZ(1)>0; END = SIZ(1)/SR; end
%    if END > START
      % sample indices, post resampling
      EIX = round(END*SR/ds);
      if SIZ(1) > 0
        EIX = min(floor(SIZ(1)/ds),EIX);
      end
      N = [SIX EIX];
%    end
%  end
  %disp(['N=',num2str(N),' MONO=',num2str(FORCEMONO),' ds=',num2str(ds)]);
  if ismp3
    [D,SR] = mp3read(FN,N,FORCEMONO,ds);
  elseif ism4a
    [D,SR] = m4aread(FN,N,FORCEMONO,ds);
  elseif iswav
    [D,SR] = wavread_downsamp(FN,N,FORCEMONO,ds);
  elseif isflac
    [D,SR] = flacread(FN,N,FORCEMONO,ds);
  else % ReadSound
    if N(2) < 0
      LEN = -1;
    else
      LEN = N(2) - N(1) + 1;
    end
    [D,SR] = ReadSound(FN,LEN,N(1)-1);
    D = D'; % returns in rows!
  end
else
  % Handlers that do not support specifying subranges of samples go here
%  if strcmp(ext, '.flac')
%    [D,SR] = flacread2(FN);
%  else
  if 1
    error(['audioread: cannot figure type of file ',FN]);
  end
  % truncate before resampling
  if START > 0 || (END > START && END < length(D)/SR)
    if END < START; END = length(D)/SR; end
    EIX = min(length(D),round(END*SR));
    D = D((round(START*SR)+1):EIX,:);
  end
end

% always mono (if not handled in mp3read etc.)
[nsamp,nchan] = size(D);
if FORCEMONO && nchan > 1
  D = mean(D,2);
end

if TARGETSR > 0 && SR ~= TARGETSR && nsamp > 1
  [p,n,e] = fileparts(FN);
%  disp(['*** audioread: resampling ',n,' from ', num2str(SR), ...
%        ' to ', num2str(TARGETSR)]);
  D = resample(D, TARGETSR, SR);
  SR = TARGETSR;
end

%disp(['audioread(',FN,',',num2str(TARGETSR),',',num2str(FORCEMONO),',',num2str(START),',',num2str(DUR),')->',num2str(size(D)),'@',num2str(SR)]);