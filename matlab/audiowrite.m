function audiowrite(D,SR,FN)
% audiowrite(D,SR,FN)
%   Write out an audio file, using wavwrite, mp3write, m4awrite, or flacwrite as
%   appropriate. 
% 2010-09-16 Dan Ellis dpwe@ee.columbia.edu

[pth,nam,ext] = fileparts(FN);

% Ensure the output directory exists
mymkdir(pth);

ext = lower(ext);
if strcmp(ext,'.mp3')
  mp3write(D,SR,FN);
elseif strcmp(ext, '.m4a') || strcmp(ext, '.aac') || strcmp(ext, '.mp4')
  m4awrite(D,SR,FN);
elseif strcmp(ext, '.flac')
  flacwrite(D,SR,FN);
elseif strcmp(ext, '.wv1') || strcmp(ext, '.wv2') || strcmp(ext, '.shn')
  % nist sphere
  % use dpwelib as a mex file
  WriteSound(FN,D',SR, 'NIST');
else
  wavwrite(D,SR,FN);
end
