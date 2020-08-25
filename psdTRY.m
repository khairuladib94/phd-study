% t=1:1:1000000;
% Fs=1;
% 
% h=spectrum.welch;
% Hx=psd(h,U,'Fs',Fs);
% 
% subplot(2,1,1)
% plot(t,U);
% subplot(2,1,2)
% plot(Hx);

% rng default
% Fs=1000;
% t=0:1/Fs:1-1/Fs;
% x=cos(2*pi*100*t)+randn(size(t));
% 
% N=length(x);
% xdft=fft(x);
% xdft=xdft(1:N/2+1);
% psdx=(1/(Fs*N))*abs(xdft).^2;
% psdx(2:end-1)=2*psdx(2:end-1);
% freq=0:Fs/length(x):Fs/2
% 
% plot(freq,10*log10(psdx))
% grid on

%remove noise/error
dnH=medfilt1(H,'omitnan');
dnD=medfilt1(D,'omitnan');
dnZ=medfilt1(Z,'omitnan');
dnF=medfilt1(F,'omitnan');

for counto=1:5
    sigH=std(dnH,'omitnan');
    muH=mean(dnH,'omitnan');
    sigD=std(dnD,'omitnan');
    muD=mean(dnD,'omitnan');
    sigZ=std(dnZ,'omitnan');
    muZ=mean(dnZ,'omitnan');
    sigF=std(dnF,'omitnan');
    muF=mean(dnF,'omitnan');
    for i=1:length(dnH)
        if dnH(i)>muH+3*sigH||dnH(i)<muH-3*sigH
            dnH(i)=NaN;
        end
        if dnD(i)>muD+3*sigD||dnD(i)<muD-3*sigD
            dnD(i)=NaN;
        end
        if dnZ(i)>muZ+3*sigZ||dnZ(i)<muZ-3*sigZ
            dnZ(i)=NaN;
        end
        if dnF(i)>muF+3*sigF||dnF(i)<muF-3*sigF
            dnF(i)=NaN;
        end
    end
end
U=abs(dnH(1:100000));
Fs=1000;
periodogram(U,rectwin(length(U)),length(U),'twosided')
