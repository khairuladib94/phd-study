% rng default
% 
% n=0:319;    %signal length is 320
% x=cos((pi/4)*n)+randn(size(n));     %Generate a sine wave signal with angular 
%                                     %frequency of pi/4 rad/sample 
%                                     %with added white noise
%                                       
% subplot(2,1,1)
% plot(x);    
% 
% pxx=pwelch(x);      %Apply Welch method. 
%                     %Using default settings: signal is divided into a
%                     %number of segments so that the length of each segment
%                     %is the longest possible but the number of segments is
%                     %not more than 8
%                     
%                     %Segments are windowed by Hamming and overlapped by 50%
%                     %by default
%                     
%                     %For this case, segment length is 71, DFT length is the
%                     %lower power of 2 than the original signal length which
%                     %is 256. This yields a frequency resolution of 2*pi/256
%                     %rad/sample
% 
% subplot(2,1,2)
% plot(10*log10(pxx));
% 
% %--------------------------------------------------------------------------
% 
% rng default
% 
% n=0:319;    %signal length is 320
% x=cos((pi/4)*n)+randn(size(n));     %Generate a sine wave signal with angular 
%                                     %frequency of pi/4 rad/sample 
%                                     %with added white noise
%                                       
% subplot(2,1,1)
% plot(x);    
% 
% pxx=pwelch(x,100,25);       %Apply Welch method. 
%                             %The second parameter is the desired window
%                             %length, means, the signal is divided into
%                             %multiple segments, each segment has 100
%                             %datapoint
%                             %The third parameter is the overlapping length
%                             %DFT length is the lower power of 2 than the 
%                             %original signal length which is 256. This 
%                             %yields a frequency resolution of 2*pi/256 rad/sample
% 
% subplot(2,1,2)
% plot(10*log10(pxx));
% 
% 
% %--------------------------------------------------------------------------
% 
% rng default
% 
% n=0:319;    %signal length is 320
% x=cos((pi/4)*n)+randn(size(n));     %Generate a sine wave signal with angular 
%                                     %frequency of pi/4 rad/sample 
%                                     %with added white noise
%                                       
% subplot(2,1,1)
% plot(x);    
% 
% pxx=pwelch(x,100,[],640);       %Apply Welch method. 
%                                 %The second parameter is the desired window
%                                 %length, means, the signal is divided into
%                                 %multiple segments, each segment has 100
%                                 %datapoint
%                                 %The third parameter is the number of
%                                 %specified DFT length
%                                 
%                                 %Since PSD estimate is one-sided, so there
%                                 %are 640/2+1 points. The highest peak lies
%                                 %on bin 81
% 
% subplot(2,1,2)
% plot(10*log10(pxx));
% xlabel('rad/sample');
% ylabel('dB');
% 
% %--------------------------------------------------------------------------

rng default

fs=1000;    %Set the sampling rate to 1kHz
t=0:1/fs:5-1/fs;    %The signal duration ranges from 0s to just before 5s
x=cos(2*pi*100*t)+randn(size(t));    %Generate the signal. The first term is a 
                                    %cosine wave (cos(2*pi*f*t)) thus
                                    %setting the signal frequency at 100 Hz
                                                                                                      
subplot(2,1,1)
plot(x);    

[pxx,f]=pwelch(x,500,300,500,fs);       %Apply Welch method. 

                                        %The second parameter is the desired window
                                        %length, 500
                               
                                        %The third parameter is the number of
                                        %overlapped samples, 300

                                        %The fourth is the DFT length, so
                                        %that 100 Hz falls direclt on DFT
                                        %bin (500 is used bcs it's the
                                        %Nyquist frequency which is fs/2?)

subplot(2,1,2)
plot(f,10*log10(pxx));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

%-----------------------------------------------------------------------------
%band-pass filtering
filtN=6;
filtfL=0.01;
filtfH=0.1;

[b,a]=butter(filtN,[filtfL filtfH]/0.5,'bandpass');

fH=filter(b,a,fg);
