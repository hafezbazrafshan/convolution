function sigConv=dirConv(sig1,sig2)
% DIRCONV  convolution using direct method
% C=DIRCONV(A,B) convolves row vectors A and B. The resulting vector is
% length LENGTH(A)+LENGTH(B)-1
% This function uses the direct method as the basis for convolution.
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

sigg=[zeros(1,flen-1) signal zeros(1,flen-1) zeros(1,flen-1)]; %zero pad signal so that it 
%doesn't run out of points while convolution
sigConv=zeros(1,clen);
filflip=fliplr(filter);
for i=flen:clen+flen-1
    sigConvv(i)=sum(sigg(i-flen+1:i).*filflip); %multiply by the flip filter 
%     and shift
end
sigConv=sigConvv(flen:end);




