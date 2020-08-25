%----------------------------------------
%Remove noise


%remove noise/error
%dnH=medfilt1(H,'omitnan');
%dnD=medfilt1(D,'omitnan');
dnZ=medfilt1(Z,'omitnan');
%dnF=medfilt1(F,'omitnan');

for i=1:5
%     sigH=std(dnH,'omitnan');
%     muH=mean(dnH,'omitnan');
%     sigD=std(dnD,'omitnan');
%     muD=mean(dnD,'omitnan');
    sigZ=std(dnZ,'omitnan');
    muZ=mean(dnZ,'omitnan');
%     sigF=std(dnF,'omitnan');
%     muF=mean(dnF,'omitnan');

    for i=1:length(dnZ)
        
%         if dnH(i)>muH+3*sigH||dnH(i)<muH-3*sigH
%             dnH(i)=NaN;
%         end
%         if dnD(i)>muD+3*sigD||dnD(i)<muD-3*sigD
%             dnD(i)=NaN;
%         end
        if dnZ(i)>muZ+3*sigZ||dnZ(i)<muZ-3*sigZ
            dnZ(i)=NaN;
        end
%         if dnF(i)>muF+3*sigF||dnF(i)<muF-3*sigF
%             dnF(i)=NaN;
%         end
        
    end
    
end

%----------------------------------------
% %FFT
Fs=1;
exp=dnZ;
nfft=length(exp);
nfft2=2^nextpow2(nfft);
ff=fft(exp,nfft2);
fff=ff(1:nfft2/2)';
xfft=Fs*(0:nfft2/2-1)/nfft2;
plot(xfft,abs(fff).*xfft);
%----------------------------------------
%Bandpass filter
% filtN=6;
% filtfL=0.01;
% filtfH=0.4;
% 
% [b,a]=butter(filtN,[filtfL filtfH]/0.5,'bandpass');
% 
% fH=filter(b,a,H(1:86400));
% fD=filter(b,a,D(1:86400));
% fZ=filter(b,a,Z(1:86400));
% fF=filter(b,a,F(1:86400)); 
% 
% subplot(2,1,1)
% plot(H);
% subplot(2,1,2)
% plot(fH);



%----------------------------------------
%FFT





%----------------------------------------
%FFT






%----------------------------------------
%FFT





%----------------------------------------
%FFT