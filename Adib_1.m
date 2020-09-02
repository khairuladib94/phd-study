%---------------Geomagnetic data analysis--------------------------------------------------------------------------
% Test 
%Data concatenation

if multyrs=='YES'     %Vertical concatenation for data of multiple years
    H=vertcat(H1,H2);
    D=vertcat(D1,D2);
    Z=vertcat(Z1,Z2);
    F=vertcat(F1,F2);
    UT1m=vertcat(UT1m1,UT1m2);
end

%Converting into daily data
dayH=reshape(H,86400,[]);
dayD=reshape(D,86400,[]);
dayZ=reshape(Z,86400,[]);
dayF=reshape(F,86400,[]);

%Nighttime data selection

tStart=0;   %Stating start and stop local time. UT will be calculated based on location
tStop=3;
if country=='PHI'
    tStart=tStart+16;
    tStop=tStop+16;
end
if country=='IND'
    tStart=tStart+17;
    tStop=tStop+17;
end    

nightH=dayH(3600*tStart+1:3600*tStop,:);
nightD=dayD(3600*tStart+1:3600*tStop,:);
nightZ=dayZ(3600*tStart+1:3600*tStop,:);
nightF=dayF(3600*tStart+1:3600*tStop,:);

%Bandpass filtering

filN=2;
filfL=0.01;
filfH=0.1;

[b,a]=butter(filN,[filfL filfH]/0.5,'bandpass');

fH=filter(b,a,nightH);
fD=filter(b,a,nightD);
fZ=filter(b,a,nightZ);
fF=filter(b,a,nightF); 


%Noise/outlier removal

dnH=medfilt1(H,'omitnan');  %Using median filter
dnD=medfilt1(D,'omitnan');
dnZ=medfilt1(Z,'omitnan');
dnF=medfilt1(F,'omitnan');

for a=1:    %Using mean+/-sigma limit
    sigH=std(dnH,'omitnan');    
    muH=mean(dnH,'omitnan');
    sigD=std(dnD,'omitnan');
    muD=mean(dnD,'omitnan');
    sigZ=std(dnZ,'omitnan');
    muZ=mean(dnZ,'omitnan');
    sigF=std(dnF,'omitnan');
    muF=mean(dnF,'omitnan');
    for b=1:length(dnH)
        if dnH(b)>muH+3*sigH||dnH(b)<muH-3*sigH
            dnH(b)=NaN;
        end
        if dnD(b)>muD+3*sigD||dnD(b)<muD-3*sigD
            dnD(b)=NaN;
        end
        if dnZ(b)>muZ+3*sigZ||dnZ(b)<muZ-3*sigZ
            dnZ(b)=NaN;
        end
        if dnF(b)>muF+3*sigF||dnF(b)<muF-3*sigF
            dnF(b)=NaN;
        end
    end
end