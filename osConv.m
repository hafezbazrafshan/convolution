function sConv=osConv(sig1,sig2,Blen,MOD)
% OSCONV  convolution using overlap and save method
% C=OSCONV(A,B,BLEN,MOD) performs convolution on A and B using overlap and
% and save method.  BLEN is the block length designated for convolution. BLEN
% should not be greater than the length of both signals. In addition BLEN
% should at least be equal to length of the shorter the signal - 1. 
% (BLEN>=FLEN-1 and FLEN=min(length(A),length(B)).
% MOD is a string that designates which method of convolution should be performed 
% in a block convolution.  MOD='direct' uses the direct method and MOD='fft' uses
% the fft.  
%The resulting vector is length LENGTH(A)+LENGTH(B)-1
%Created by Mohammadhafez Bazrafshan
if nargin<2
    error('Signals Not defined');
end
if nargin<3
    error('Block Length Not defined');
end

if nargin<4
    error('Mode of conv not defined, ''direct'' for direct convolution and ''fft'' for fast');
end

[A AA]=size(sig1);
[B BB]=size(sig2);
log1=(A~=1) & (AA~=1);
log2=(B~=1) & (BB~=1);
if  log1==1
    error('signal 1 should be one dimensional');
end
if  log2==1
    error('signal 2 should be one dimensional');
  
end
        
if length(sig1)>length(sig2)
    filter=sig2;
    signal=sig1;
else
    filter=sig1;
    signal=sig2;
end


flen=length(filter);
slen=length(signal);
clen=flen+slen-1;  %length of convolution


if Blen>slen
    error('Block length should be less than signal length');
end
if Blen==0
    error('Block length Should not be zero');
end
if Blen<flen-1
    error('Block length should be at least filter length-1');
end




r=rem(slen,Blen);
if r~=0
sigg=[signal zeros(1,Blen-r)];
else 
    sigg=signal;
end


sslen=length(sigg);
Nbins=sslen./Blen;



switch MOD
case 'direct'
sConvv=zeros(1,Nbins*Blen+flen-1);
sigg=[zeros(1,flen-1) sigg];
idx1=1;

for i=1:Nbins
    
    ss2=sigg(idx1:idx1+Blen-1+flen-1);  %taking out one bin = N1 points
    cnv=dirConv(ss2,filter);
    sConvv(idx1:idx1+Blen-1+flen-1)=cnv(flen:end);
    idx1=idx1+Blen;
end


    
sConv=sConvv(1:clen); 


case 'fft'

sConvv=zeros(1,Nbins*Blen+flen-1);
sigg=[zeros(1,flen-1) sigg];
idx1=1;

for i=1:Nbins
    
    ss2=sigg(idx1:idx1+Blen-1+flen-1);  %taking out one bin = N1 points
    cnv=fftConv(ss2,filter);
    sConvv(idx1:idx1+Blen-1+flen-1)=cnv(flen:end);
    idx1=idx1+Blen;
end

    
sConv=sConvv(1:clen); 
end

    
    
    
    


