function fConv=fftConv(sig1,sig2)
% FFTCONV  convolution using fast fourier transform
% C=FFTCONV(A,B) convolves row vectors A and B. The resulting vector is
% length LENGTH(A)+LENGTH(B)-1
% This function uses the FFT as the basis for convolution.
% Created by Mohammadhafez Bazrafshan
if nargin<2
    error('Input Arguments Not defined');
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

fConv=real(ifft(fft(signal,clen).*fft(filter,clen)));