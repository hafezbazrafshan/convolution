%EE-5163, Digital Signal Processing, ART GRIGORYAN, UTSA-2013
%Project: Linear Filtering of Signals
%Student name: Mohammadhafez Bazrafshan
%Banner ID: @01347682
%Due March 7th, 2013. Extended to March 19th, 2013



clear all;
clear;
clc;
close all;



%***********************PART 1 reading the signal*************************%
[fn,Fs,nbits]=wavread('mike.wav');
fn=fn'; %column vector to row vector (I'm more comfortable);
% Fs - sampling rate  22050 Hz
N=length(fn);       % 854128 
N_insec = N/Fs;    % 38.7360 sec
% sound(X,Fs);
FN=abs(fft(fn));
% ************************************************************************%


% ************************PART 2 NORMAL DISTRIBUTION**********************%
std=0.01; %standard deviation of noise
meaN=0;   %mean of Noise
nN=meaN+std*randn(1,N);   %generate random noise with normal distribution
gn=fn+nN;  %add noise to the original signal
GN=abs(fft(gn));
% ************************************************************************%
B=5;  %5 seconds partitioning
or1=1; %position of zero (origin) in gn 

figure(1)
subplot(2,1,1);
plot(linspace(0,N_insec,N),fn)
grid on
xlabel('time (seconds)');
ylabel('Amplitude');
title('speech signal');
axis([0 N_insec -0.5 0.5]);
subplot(2,1,2);
plot(linspace(0,N_insec,N),gn);
grid on;
xlabel('time (seconds)');
ylabel('Amplitude');
title('Noise added to speech signal');
axis([0 N_insec -0.5 0.5]);

figure(2)
subplot(2,1,1);
fvec=linspace(-Fs/2,Fs/2,N)./1000;
plot(fvec,fftshift(FN));
xlabel('Frequency (kHz)');
ylabel('FFT of speech signal');

subplot(2,1,2);
fvec=linspace(-Fs/2,Fs/2,N)./1000;
plot(fvec,fftshift(GN));
xlabel('Frequency (kHz)');
ylabel('FFT of Noisy signal');

figure(3)
subplot(2,1,1);
fvec=[-0.1:Fs./N:0.1];
plot(fvec,fftshift(FN(1:length(fvec))));
xlabel('Frequency (Hz)');
ylabel('FFT of speech signal');
title('Zoomed In');
subplot(2,1,2);
fvec=[-0.1:Fs./N:0.1];
plot(fvec,fftshift(GN(1:length(fvec))));
xlabel('Frequency (Hz)');
ylabel('FFT of Noisy signal');

% ************************PART 3 TRIANGLE FILTER**************************%
hn=[1 1 2 2 3 3 4 4 5 5 6 6 8 7 7 6 6 5 5 4 4 3 3 2 2 1 1]; %triangle
NH=length(hn); % length of the filter
or2=13;  %position of zero(origin) in hn
NH_insec=NH./Fs;
figure(4);
plot([-12:14],hn,'-');
xlabel('time index (nT_s)');
ylabel('amplitude');
title('Triangle Filter');
grid on;

figure(5);
plot(linspace(-Fs./2,Fs./2,NH)./1000,abs(fftshift(fft(hn))));
title('FFT of Triangle Filter');
xlabel('Frequency (kHz)');
ylabel('amplitude');
grid on;
hn=hn./sum(hn); %normalizing hn; if we want to listen to sound
% ************************************************************************%

% ***********************PART 4 Direct Convolution************************%
yn4=zeros(1,length(fn)+length(hn)-1);
yn4=dirConv(fn,hn);
tic
ygn4=dirConv(gn,hn);
toc
%Elapsed time is 4.476943 seconds.
figure(6);
subplot(2,1,1);
plot(linspace(-NH_insec/2,N_insec+NH_insec/2,length(yn4)),yn4);
axis([-NH_insec N_insec+NH_insec/2 -0.5 0.5]);
xlabel('time (seconds)');
ylabel('Amplitude');
title('Direct Convolution of Speech Signal');
grid on;
subplot(2,1,2);
plot(linspace(-NH_insec/2,N_insec+NH_insec/2,length(ygn4)),ygn4);
axis([-NH_insec N_insec+NH_insec/2 -0.5 0.5]);
xlabel('time (seconds)');
ylabel('Amplitude');
title('Direct Convolution of Noisy Signal');
grid on;
dirErr=sqrt(sum((ygn4(or2:or2+N-1)-fn).^2))./N;


%***********************PART 5 FFT Based Convolution**********************%
yn5=fftConv(fn,hn);
tic
ygn5=fftConv(gn,hn);
toc
%Elapsed time is 0.323113 seconds.
figure(7);
subplot(2,1,1);
plot(linspace(-NH_insec/2,N_insec+NH_insec/2,length(yn5)),yn5);
axis([-NH_insec N_insec+NH_insec/2 -0.5 0.5]);
xlabel('time (seconds)');
ylabel('Amplitude');
title('FFT Based Convolution of Speech Signal');
grid on;
subplot(2,1,2);
plot(linspace(-NH_insec/2,N_insec+NH_insec/2,length(ygn5)),ygn5);
axis([-NH_insec N_insec+NH_insec/2 -0.5 0.5]);
xlabel('time (seconds)');
ylabel('Amplitude');
title('FFT Based Convolution of Noisy Signal');
grid on;
fftErr=sqrt(sum((ygn5(or2:or2+N-1)-fn).^2))./N;


%***********************PART 6.1 Overlap-Save Convolution*****************%
blen=5*Fs;  %block length
yn61=osConv(fn,hn,blen,'direct');
tic
ygn61=osConv(gn,hn,blen,'direct');
toc
%elapsed time;
figure(8);
subplot(2,1,1);
plot(linspace(-NH_insec/2,N_insec+NH_insec/2,length(yn61)),yn61);
axis([-NH_insec N_insec+NH_insec/2 -0.5 0.5]);
xlabel('time (seconds)');
ylabel('Amplitude');
title('Overlap-Save Direct Convolution of Speech Signal');
grid on;
subplot(2,1,2);
plot(linspace(-NH_insec/2,N_insec+NH_insec/2,length(ygn61)),ygn61);
axis([-NH_insec N_insec+NH_insec/2 -0.5 0.5]);
xlabel('time (seconds)');
ylabel('Amplitude');
title('Overlap-Save Direct Convolution of Noisy Signal');
grid on;
osErr=sqrt(sum((ygn61(or2:or2+N-1)-fn).^2))./N;






%***********************PART 6.2 Overlap-Save Convolution*****************%
blen=5*Fs;  %block length
yn62=oaConv(fn,hn,blen,'direct');
tic
ygn62=oaConv(gn,hn,blen,'direct');
toc
%elapsed time;
figure(9);
subplot(2,1,1);
plot(linspace(-NH_insec/2,N_insec+NH_insec/2,length(yn62)),yn62);
axis([-NH_insec N_insec+NH_insec/2 -0.5 0.5]);
xlabel('time (seconds)');
ylabel('Amplitude');
title('Overlap-Add Direct Convolution of Speech Signal');
grid on;
subplot(2,1,2);
plot(linspace(-NH_insec/2,N_insec+NH_insec/2,length(ygn62)),ygn62);
axis([-NH_insec N_insec+NH_insec/2 -0.5 0.5]);
xlabel('time (seconds)');
ylabel('Amplitude');
title('Overlap-Add Direct Convolution of Noisy Signal');
grid on;
oaErr=sqrt(sum((ygn62(or2:or2+N-1)-fn).^2))./N;



%************PART 6.3 Overlap-Save and Overlap-Add FFT Convolution********%
yn63a=osConv(fn,hn,blen,'fft');
tic
ygn63a=osConv(gn,hn,blen,'fft');
toc
%elapsed time;
figure(10);
subplot(2,1,1);
plot(linspace(-NH_insec/2,N_insec+NH_insec/2,length(yn63a)),yn63a);
axis([-NH_insec N_insec+NH_insec/2 -0.5 0.5]);
xlabel('time (seconds)');
ylabel('Amplitude');
title('Overlap-Save FFT Convolution of Speech Signal');
grid on;
subplot(2,1,2);
plot(linspace(-NH_insec/2,N_insec+NH_insec/2,length(ygn63a)),ygn63a);
axis([-NH_insec N_insec+NH_insec/2 -0.5 0.5]);
xlabel('time (seconds)');
ylabel('Amplitude');
title('Overlap-Save FFT Convolution of Noisy Signal');
grid on;
osfftErr=sqrt(sum((ygn63a(or2:or2+N-1)-fn).^2))./N;


yn63b=oaConv(fn,hn,blen,'fft');
tic
ygn63b=oaConv(gn,hn,blen,'fft');
toc
%elapsed time;
figure(11);
subplot(2,1,1);
plot(linspace(-NH_insec/2,N_insec+NH_insec/2,length(yn63b)),yn63b);
axis([-NH_insec N_insec+NH_insec/2 -0.5 0.5]);
xlabel('time (seconds)');
ylabel('Amplitude');
title('Overlap-Add FFT Convolution of Speech Signal');
grid on;
subplot(2,1,2);
plot(linspace(-NH_insec/2,N_insec+NH_insec/2,length(ygn63b)),ygn63b);
axis([-NH_insec N_insec+NH_insec/2 -0.5 0.5]);
xlabel('time (seconds)');
ylabel('Amplitude');
title('Overlap-Add FFT Convolution of Noisy Signal');
grid on;
oafftErr=sqrt(sum((ygn63b(or2:or2+N-1)-fn).^2))./N;


% 
% Elapsed time is 4.618433 seconds.
% Elapsed time is 0.346697 seconds.
% Elapsed time is 4.722118 seconds.
% Elapsed time is 5.087066 seconds.
% Elapsed time is 0.293950 seconds.
% Elapsed time is 0.376626 seconds.