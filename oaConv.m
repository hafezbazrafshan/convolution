function aConv=oaConv(sig1,sig2,Blen,MOD)
% OACONV  convolution using overlap and add method
% C=OACONV(A,B,BLEN,MOD) performs convolution on A and B using overlap and
% and add method.  BLEN is the block length designated for convolution. BLEN
% should not be greater than the length of both signals. 
% MOD is a string that designates which method of convolution should be performed 
% in a block convolution.  MOD='direct' uses the direct method and MOD='fft' uses
% the fft.  
%The resulting vector is length LENGTH(A)+LENGTH(B)-1
%Created by Mohammadhafez Bazrafshan

if nargin<2
    error('Input Arguments Not defined');
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




r=rem(slen,Blen);
if r~=0
sigg=[signal zeros(1,Blen-r)];
else 
    sigg=signal;
end


switch MOD
    case 'direct'
sslen=length(sigg);
Nbins=sslen./Blen;
aConvv=zeros(1,Nbins*Blen+flen-1);

for i=1:Nbins
    idx1=(i-1)*Blen+1;
    idx2=i*Blen;
    Sig=sigg(idx1:idx2);
    aConvv(idx1:idx2+flen-1)=aConvv(idx1:idx2+flen-1)+dirConv(Sig,filter);
end
aConv=aConvv(1:clen);

    case 'fft'
        sslen=length(sigg);
Nbins=sslen./Blen;
aConvv=zeros(1,Nbins*Blen+flen-1);

for i=1:Nbins
    idx1=(i-1)*Blen+1;
    idx2=i*Blen;
        Sig=sigg(idx1:idx2);
    aConvv(idx1:idx2+flen-1)=aConvv(idx1:idx2+flen-1)+fftConv(Sig,filter);
end
aConv=aConvv(1:clen);

end
        
