function [lev,fso]=activlev(sp,fs,mode)
%ACTIVLEV Measure active speech level as in ITU-T P.56 [LEV,FSO]=(SP,FS,MODE)
%
%Inputs: SP is the speech signal (with better than 20dB SNR)
%        FS is the sample frequency in Hz (see also FSO below)
%        MODE is any combination of the following:
%            r - raw: do not apply the input filter 0.2 to 5.5 kHz
%            a - all: use all the samples rather than subsampling to 694Hz
%            d - give outputs in dB rather than power
%            l - give the "long term level" as well as the active speech level
%Outputs:
%        LEV gives the speech level in units of power (or dB if mode='d')
%            if mode='l' is specified, LEV is a row vector with the "long term
%            level" as its second element (this is just the mean power)
%        FSO is a column vector of intermediate information that allows
%            you to process a speech signal in chunks. Thus:
%
%            fso=fs; for i=1:inc:nsamp, [lev,fso]=activlev(sp(i:i+inc-1),fso,mode); end
%
%            is equivalent to:                lev=activlev(sp(1:nsamp),fs,mode)
%
%            but is much slower. The two methods will not give identical results
%            because they will use slightly different threshods.

%For completeness we list here the contents of the FSO array:
%
%     1 : sample frequency
%     2 : subsampling increment
%     3 : hangover increment
%   4-6 : smoothing filter coefs
%  7-12 : 200Hz HP filter numerator
% 13-18 : 200Hz HP filter denominator
% 19-24 : 5.5kHz LP filter numerator
% 25-30 : 5.5kHz LP filter denominator
% 31-33 : count,sum and sum of squares of filter speech
%    34 : upper limit of largest amplitude bin
%    35 : samples to skip for subsampling
% 36-37 : smoothing filter state
% 38-42 : 200Hz HP filter state
% 43-47 : 5.5kHz LP filter state
% 48-55 : thresholded sample counts
% 56-63 : hangover counts

%      Copyright (C) Mike Brookes 1998
%
%      Last modified Thu Apr 30 17:28:29 1998
%
%   VOICEBOX home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   ftp://prep.ai.mit.edu/pub/gnu/COPYING-2.0 or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3 mode='0'; else mode = [mode '0']; end
lev=zeros(1,1+any(mode=='l'));
thf=0.5;
nb=8;
fso=[fs;zeros(64-1-length(fs),1)];
if length(fs)==1
   ni=max(floor(fs*all(mode~='a')/694),1);
   ti=ni/fs;
   nh=ceil(0.2/ti)+1;
   g=exp(-ti/0.03);
   % s-plane zeros and poles of high pass 5'th order chebychev2 filter with -0.25dB at w=1
   szp=[0.37843443673309i 0.23388534441447i; -0.20640255179496+0.73942185906851i -0.54036889596392+0.45698784092898i]; 
   szp=[[0; -0.66793268833792] szp conj(szp)];
   zl=2./(1-szp*tan(200*pi/fs))-1;
   al=real(poly(zl(2,:)));
   bl=real(poly(zl(1,:)));
   sw=1-2*rem(0:5,2).';
   bl=bl*(al*sw)/(bl*sw);
   fso(2:19-1)=[ni;nh;[1; -2*g; g^2]/(1-2*g+g^2);bl(:);al(:)];
   if fs>14000
      zh=2./(szp/tan(5500*pi/fs)-1)+1;
      ah=real(poly(zh(2,:)));
      bh=real(poly(zh(1,:)));
      bh=bh*sum(ah)/sum(bh);
      fso(19:19+11)=[bh(:);ah(:)];
   end
end
% need to calculate offset correctly
% subsample the signal
ns=length(sp);
if ns
   % input filter goes here
   if all(mode~='r')
      [sp,fso(38:42)]=filter(fso(7:12),fso(13:18),sp(:),fso(38:42));
      if fso(25)
         [sp,fso(43:47)]=filter(fso(19:24),fso(25:30),sp,fso(43:47));
      end
   end
   sp=abs(sp(1+fso(35):fso(2):ns));
   na=length(sp);
   fso(35)=fso(35)+na*fso(2)-ns;
   if (na)
      fso(31:33)=fso(31:33)+[na;sum(sp);sp.'*sp];
      [sp,fso(36:37)]=filter(1,fso(4:6),sp,fso(36:37));
      % determine new peak
      th=max(sp);
      s=zeros(nb,1);
      h=s;
      if fso(34)>0
         if th>fso(34)
            nt=floor(log(th/fso(34))/log(thf));
            th=fso(34)*thf^nt;
            if nt>-nb
               s(1-nt:nb)=fso(48:47+nb+nt);
               h(1-nt:nb)=fso(56:55+nb+nt);
            end
         else
            th=fso(34);
            s=fso(48:55);
            h=fso(56:63);
         end
      end
      fso(34)=th;
      nh=fso(3);
      if nh<na
         nm=na-nh+2;
         for i=1:nb
            th=th*thf;
            b=sp>th;
            np=h(i);
            s(i)=s(i)+np+sum(cumsum([sum(b(1:np+1));b(2+np:na)]-[zeros(nh-np,1);b(1:na-nh)])>0);
            ff=find([1;b(nm:na)]);
            nm=nm-1+ff(end);
            h(i)=nm+nh-na-2;;
         end
      else
         for i=1:nb
            th=th*thf;
            b=sp>th;
            np=h(i);
            ff=find(b);
            if length(ff)
               h(i)=ff(end)+nh-na-1;
               s(i)=s(i)+na+min(0,np+1-ff(1));
            else
               h(i)=max(0,np-na);
               s(i)=s(i)-h(i)+np;
            end
         end
      end
      fso(48:63)=[s;h];
   end
end
ssq=fso(33);
if ssq>0
   s=fso(48:55);
   thresh=(fso(34)*thf.^(1:nb).').^2;
   % 15.9 dB = 38.9 power ratio
   mar=38.9;
   st=s.*thresh*mar/ssq;
   ff=find(st>1);
   if length(ff)
      nf=ff(end);
      ls=log(s(nf)/s(min(nb,nf+1)));
      lev(1)=ssq/s(nf)*st(nf)^(ls/(ls-2*log(thf)));
   else
      [ls,nf]=max(st);
      lev(1)=ssq/s(nf);
   end
   if length(lev)>1 lev(2)=ssq/fso(31); end
end
if any(mode=='d')
   lev=10*log10(max(lev,1e-200));
end
