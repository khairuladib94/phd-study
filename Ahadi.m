%Earthquake----------------------------------------------------------------

%Seismic Index

% if Mw<6.2                   %Converting Mw (moment magnitude) to Ms (surface wave magnitude)
%     Ms=(Mw-2.07)/0.67;
% end
% 
% if Mw>~6.2
%     Ms=(Mw-0.08)/0.99;
% end
% 
% Ks=(1+R^(-Ms/2))^(-2.33)+(10^(0.75*Ms)/10*R);            %Calculating Seismic index

%Geomagnetic data----------------------------------------------------------

%Remove noise/outlier

dnH=medfilt1(H,'omitnan');
dnD=medfilt1(D,'omitnan');
dnZ=medfilt1(Z,'omitnan');

for i=1:5
    
    sigH=std(dnH,'omitnan');
    muH=mean(dnH,'omitnan');
    sigD=std(dnD,'omitnan');
    muD=mean(dnD,'omitnan');
    sigZ=std(dnZ,'omitnan');
    muZ=mean(dnZ,'omitnan');
    
    for i=1:length(dnZ)
        
        if dnH(i)>muH+3*sigH||dnH(i)<muH-3*sigH
            dnH(i)=NaN;
        end
        if dnD(i)>muD+3*sigD||dnD(i)<muD-3*sigD
            dnD(i)=NaN;
        end
        if dnZ(i)>muZ+3*sigZ||dnZ(i)<muZ-3*sigZ
            dnZ(i)=NaN;
        end
        
    end
    
end

% figure(1)
% subplot(3,1,1)
% plot(dnH)
% subplot(3,1,2)
% plot(dnD)
% subplot(3,1,3)
% plot(dnZ)
% 
% figure(2)
% subplot(3,1,1)
% plot(dtH)
% subplot(3,1,2)
% plot(dtD)
% subplot(3,1,3)
% plot(dtZ)


%Daily data

dayH=reshape(dnH,86400,[]);
dayD=reshape(dnD,86400,[]);
dayZ=reshape(dnZ,86400,[]);

%Nighttime data 2200 - 0300 hrs LT. Insert time in UT

tStart=15;      %Indonesia timezone -> UT
tStop=20;

nightH=dayH(3600*tStart:3600*tStop-1,:);
nightD=dayD(3600*tStart:3600*tStop-1,:);
nightZ=dayZ(3600*tStart:3600*tStop-1,:);

%Remerging data
MergeH=reshape(nightH,[],1);
MergeD=reshape(nightD,[],1);
MergeZ=reshape(nightZ,[],1);

% %FFT
% Fs=1;
% % nfft=length(MergeZ(1:102);
% nfft2=2^22;
% window=hamming(nfft2);
% ff=fft(MergeZ(1:2^22),nfft2);
% fff=ff(1:nfft2/2)';
% xfft=Fs*(0:nfft2/2-1)/nfft2;
% plot(xfft,abs(fff).*xfft);


fs=1; 
t=0:1/fs:length(Hsel);

                                                                                                      
subplot(4,1,1)
plot(fg);    

dtZ=detrend(fg)
subplot(4,1,2)
plot(dtZ);

[pxx,f]=pwelch(H_filt,1024,[],[],fs);                %Apply Welch method. 

                                                %The second parameter is the desired window
                                                %length, 500

                                                %The third parameter is the number of
                                                %overlapped samples, 300

                                                %The fourth is the DFT length, so
                                                %that 100 Hz falls direclt on DFT
                                                %bin (500 is used bcs it's the
                                                %Nyquist frequency which is fs/2?)

% subplot(4,1,3)
figure(4)
plot(f,10*log10(pxx));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');



