%Input
stn='CDO';                              %abbreviation of the station name, 3 letters. e.g: 'DAV'
station='Cagayan De Oro';                        %full name of the station name. e.g: Davao
country='PHI';                          %abbreviation of country name, 3 letters. e.g: 'IND','PHI'
minmag=5.0;                             %minimum of earthquake magnitude. e.g: 5.0
maxepidis=300;                          %maximum epicentral distance (in km). e.g 150. Insert 'suggest' to get suggestion for maximum epicentral distance as described by Uyeda, 2003

%----------------------------------------------------------------------------------------------------------------------------------------------
%Earthquake occurrances

year=min(datevec(UT1m));
year=year(1,1);
if maxepidis=='suggest'
    maxepidis=(minmag-4.5)/0.025;
end
%check station
if stn=='TGG'
    stnnum=1;
end
if stn=='MUT'
    stnnum=2;
end
if stn=='LGZ'
    stnnum=3;
end
if stn=='CEB'
    stnnum=4;
end
if stn=='CDO'
    stnnum=5;
end
if stn=='DAV'
    stnnum=6;
end
if stn=='GSI'
    stnnum=7;
end
if stn=='SCN'
    stnnum=8;
end
if stn=='LWA'
    stnnum=9;
end
if stn=='PTN'
    stnnum=10;
end
if stn=='MND';
stnnum=11;
end
if stn=='BIK'
    stnnum=12;
end
if stn=='JYP'
    stnnum=13;
end
if stn=='PRP'
    stnnum=14;
end
if stn=='KPG'
    stnnum=15;
end
    
latlonstn=[StnIndPhillip(stnnum,2) StnIndPhillip(stnnum,3)];
    
%check leap year
if rem(year,4)==0
    days=1:1:366;
    numdays=366;
else
    days=1:1:365;
    numdays=365;
end

%check earthquakes
if country=='IND'
    eqvar=EqIndo20112017;
    else if country=='PHI'
    eqvar=EqPhillip20112017;
    end
end

toteq=length(eqvar);
dataeq=NaN(100,4);
dataeqnum=1;
for numeq=1:toteq
    if eqvar(numeq,1)==year
        if  eqvar(numeq,10)>=minmag
            latloneq=[eqvar(numeq,7) eqvar(numeq,8)];
            epidis=deg2km(distance('gc',latlonstn,latloneq));
            dateeq=datenum([eqvar(numeq,1) eqvar(numeq,2) eqvar(numeq,3)]);
            mag=eqvar(numeq,10);
            daynum=day(datetime([eqvar(numeq,1) eqvar(numeq,2) eqvar(numeq,3)]),'dayofyear');
 
            if epidis<=maxepidis
                dataeq(dataeqnum,1)=dateeq;
                dataeq(dataeqnum,2)=daynum;
                dataeq(dataeqnum,3)=epidis;
                dataeq(dataeqnum,4)=mag;
                dataeqnum=dataeqnum+1;
            end
        end
    end
end

dataeq(any(isnan(dataeq),2),:)=[];
if size(dataeq,1)>1
    dataeq=flip(dataeq);
end

for countj=1:3
    for countg=1:size(dataeq,1)
        if countg~=1
         if dataeq(countg,2)==dataeq(countg-1,2)
             if dataeq(countg,4)<dataeq(countg-1,4)
                dataeq(countg,:)=NaN;
             else
                dataeq(countg-1,:)=NaN;
             end
         end
        end
    end
    dataeq(any(isnan(dataeq),2),:)=[];
end

eqstring=string(zeros(size(dataeq,1),1));
eqletter={'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','A1','B1','C1','D1','E1','F1','G1','H1','I1','J1','K1','L1','M1','N1','O1','P1','Q1','R1','S1','T1','U1','V1','W1','X1','Y1','Z1','A2','B2','C2','D2','E2','F2','G2','H2','I2','J2','K2','L2','M2','N2','O2','P2','Q2','R2','S2','T2','U2','V2','W2','X2','Y2','Z2'};
eql=zeros(size(dataeq,1),1);
eql=(string(eqletter(1:(size(dataeq,1)))))';

eqlbl=eql;
for counth=1:size(dataeq,1)
    if counth~=1
     if (dataeq(counth,2)-dataeq(counth-1,2))<=2
        eqlbl(counth)=sprintf('%s%s',eqlbl(counth-1),eqlbl(counth));
        eqlbl(counth-1)=' ';
     end
    end
end

for countf=1:size(dataeq,1)
   date=datetime(datevec(dataeq(countf,1)),'Format','dd/MM');
   eqstring(countf,1)=sprintf('%s, M=%.1f, %.0f km',date,dataeq(countf,4),dataeq(countf,3));
end

if numel(eql)>0
    eqstring=strcat('(',eql,')',eqstring);
end

%----------------------------------------------------------------------------------------------------------------------------------------------
%Earthquakes Index

Eq20112017=vertcat(EqIndo20112017,EqPhillip20112017);
Eq20112017=sortrows(Eq20112017,3);
Eq20112017=sortrows(Eq20112017,2);
Eq20112017=sortrows(Eq20112017,1);

EqInd=zeros(numdays,1);
countx=1;

for county=1:length(Eq20112017)
    
    if Eq20112017(county,1)==year
        
        daynum=day(datetime([Eq20112017(county,1) Eq20112017(county,2) Eq20112017(county,3)]),'dayofyear');
        latloneq=[Eq20112017(county,7) Eq20112017(county,8)];
        dis=deg2km(distance('gc',latlonstn,latloneq));
        mag=Eq20112017(county,10);
        idx=mag^3/dis^2;
        
        if (dis<=1000) && (daynum==countx)
            EqInd(countx)=EqInd(countx)+idx;
        end
        
        if (dis<=1000) && (daynum~=countx)
            countx=daynum;
            EqInd(countx)=EqInd(countx)+idx;
        end
        
    end
end
    
%----------------------------------------------------------------------------------------------------------------------------------------------
%Geomagnetic field

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

%reshape arrays into daily matrices
dayH=reshape(dnH,86400,[]);
dayD=reshape(dnD,86400,[]);
dayZ=reshape(dnZ,86400,[]);
dayF=reshape(dnF,86400,[]);

%nighttime data only 0000 - 0400 hrs LT. Insert time in UT
tStart=16;
tStop=20;
nightH=dayH(3600*tStart+1:3600*tStop,:);
nightD=dayD(3600*tStart+1:3600*tStop,:);
nightZ=dayZ(3600*tStart+1:3600*tStop,:);
nightF=dayF(3600*tStart+1:3600*tStop,:);

%band-pass filtering
filtN=6;
filtfL=0.01;
filtfH=0.1;

[b,a]=butter(filtN,[filtfL filtfH]/0.5,'bandpass');

fH=filter(b,a,nightH);
fD=filter(b,a,nightD);
fZ=filter(b,a,nightZ);
fF=filter(b,a,nightF); 

%power spectral densities
deltaf=filtfH-filtfL;
sH=((abs(fH)).^2)/(2*pi*deltaf);
sD=((abs(fD)).^2)/(2*pi*deltaf);
sZ=((abs(fZ)).^2)/(2*pi*deltaf);
sF=((abs(fF)).^2)/(2*pi*deltaf);

%reshape arrays into monthly matrices
smonH=reshape(sH,[],12);
smonD=reshape(sD,[],12);
smonZ=reshape(sZ,[],12);
smonF=reshape(sF,[],12);

%daily mean
daymeanH=mean(sH,1,'omitnan');
daymeanD=mean(sD,1,'omitnan');
daymeanZ=mean(sZ,1,'omitnan');
daymeanF=mean(sF,1,'omitnan');

%monthly mean
monmeanH=mean(smonH,1,'omitnan');
monmeanD=mean(smonD,1,'omitnan');
monmeanZ=mean(smonZ,1,'omitnan');
monmeanF=mean(smonF,1,'omitnan');

%monthly standard deviation
monstdH=std(smonH,0,1,'omitnan');
monstdD=std(smonD,0,1,'omitnan');
monstdZ=std(smonZ,0,1,'omitnan');
monstdF=std(smonF,0,1,'omitnan');

%normalisation

counta=1;
monI=1;
for dayI=1:numdays
        
       normH(1,dayI)=(daymeanH(1,dayI)-monmeanH(1,monI))/monstdH(1,monI);
       normD(1,dayI)=(daymeanD(1,dayI)-monmeanD(1,monI))/monstdD(1,monI);
       normZ(1,dayI)=(daymeanZ(1,dayI)-monmeanZ(1,monI))/monstdZ(1,monI);
       normF(1,dayI)=(daymeanF(1,dayI)-monmeanF(1,monI))/monstdF(1,monI);
       counta=counta+1;
       
       if counta==32
        monI=monI+1;
        counta=1;
       end
  
end

%Z/H and Z/D polarisation ratio
ZH=normZ./normH;
ZD=normZ./normD;

normG=sqrt(normH.^2+normD.^2);
ZG=normZ./normG;


%mu+/-2sigma
totmeanZH=mean(ZH,'omitnan');
totmeanZD=mean(ZD,'omitnan');
totstdZH=std(ZH,'omitnan');
totstdZD=std(ZD,'omitnan');
mu2psigZH=ones(size(days))*(totmeanZH+2*totstdZH);
mu2msigZH=ones(size(days))*(totmeanZH-2*totstdZH);
mu2psigZD=ones(size(days))*(totmeanZD+2*totstdZD);
mu2msigZD=ones(size(days))*(totmeanZD-2*totstdZD);

totmeanZG=mean(ZG,'omitnan');
totstdZG=std(ZG,'omitnan');
mu2psigZG=ones(size(days))*(totmeanZG+2*totstdZG);
mu2msigZG=ones(size(days))*(totmeanZG-2*totstdZG);

%----------------------------------------------------------------------------------------------------------------------------------------------
%Plotting

subplot(4,1,1);
titlestr = sprintf('%s station, %d | Earthquakes: Mw_m_i_n=%.1f, d_m_a_x=%.0f km',station,year,minmag,maxepidis);
plot(days,ZH,'r',days,mu2psigZH,'k--',days,mu2msigZH,'k--');
hold on
plot(days,movmean(ZH,10),'b');
hold off
axis([1 numdays -inf inf]);
set(gca,'XTick',[dataeq(:,2)],'XTickLabel',{eqlbl});
if (min(mu2msigZH)<0 && min(mu2psigZH)>0)
    set(gca,'YTick',[min(mu2msigZH) 0 min(mu2psigZH)],'YTickLabel',{'\mu-2\sigma','0','\mu+2\sigma'});
    else if (min(mu2msigZH)<0 && min(mu2psigZH)<0)
        set(gca,'YTick',[min(mu2msigZH) min(mu2psigZH) 0],'YTickLabel',{'\mu-2\sigma','\mu+2\sigma','0'});
        else set(gca,'YTick',[0 min(mu2msigZH) min(mu2psigZH)],'YTickLabel',{'0','\mu-2\sigma','\mu+2\sigma'});
    end
end   
ylabel('^{Z}/_{H}');
title(titlestr);

subplot(4,1,2);
plot(days,ZD,'r',days,mu2psigZD,'k--',days,mu2msigZD,'k--');
hold on
plot(days,movmean(ZD,10),'b');
hold off
axis([1 numdays -inf inf]);
set(gca,'XTick',[dataeq(:,2)],'XTickLabel',{eqlbl});
if (min(mu2msigZD)<0 && min(mu2psigZD)>0)
    set(gca,'YTick',[min(mu2msigZD) 0 min(mu2psigZD)],'YTickLabel',{'\mu-2\sigma','0','\mu+2\sigma'});
    else if (min(mu2msigZD)<0 && min(mu2psigZD)<0)
        set(gca,'YTick',[min(mu2msigZD) min(mu2psigZD) 0],'YTickLabel',{'\mu-2\sigma','\mu+2\sigma','0'});
        else set(gca,'YTick',[0 min(mu2msigZD) min(mu2psigZD)],'YTickLabel',{'0','\mu-2\sigma','\mu+2\sigma'});
    end
end        
ylabel('^{Z}/_{D}');

% subplot(4,1,3);
% plot(days,ZG,'r',days,mu2psigZG,'k--',days,mu2msigZG,'k--');
% hold on
% plot(days,movmean(ZG,10),'b');
% hold off
% axis([0 numdays -inf inf]);
% set(gca,'XTick',[dataeq(:,2)],'XTickLabel',{eqlbl});
% if (min(mu2msigZG)<0 && min(mu2psigZG)>0)
%     set(gca,'YTick',[min(mu2msigZG) 0 min(mu2psigZG)],'YTickLabel',{'\mu-2\sigma','0','\mu+2\sigma'});
%     else if (min(mu2msigZG)<0 && min(mu2psigZG)<0)
%         set(gca,'YTick',[min(mu2msigZG) min(mu2psigZG) 0],'YTickLabel',{'\mu-2\sigma','\mu+2\sigma','0'});
%         else set(gca,'YTick',[0 min(mu2msigZG) min(mu2psigZG)],'YTickLabel',{'0','\mu-2\sigma','\mu+2\sigma'});
%     end
% end        
% ylabel('^{Z}/_{G}');

subplot(4,1,3)
plot(days,EqInd);
axis([1 numdays -inf inf]);
set(gca,'XTick',[dataeq(:,2)],'XTickLabel',{eqlbl});
ylabel('Earthquake Index');

Kp=zeros(numdays,1);
Dst=zeros(numdays,1);
counte=1;
for countd=1:length(Ind20112017)    
    if Ind20112017(countd,1)==year
         Kp(counte,1)=Ind20112017(countd,5);
         Dst(counte,1)=Ind20112017(countd,6);
         counte=counte+1;
    end  
end

subplot(4,1,4);
bar(days,Kp,'blue');
hold on
bar(days,Dst,'red');
hold off
xlabel('Days');
ylabel('Dst             \SigmaK_{p}');
legend('\SigmaK_{p}','Dst','location','southeast');
axis([1 numdays -inf inf]);

Details=eqstring;
dim = [0.905 0.725 0.2 0.2];
if numel(Details)~=0
    annotation('textbox',dim,'String',Details,'FitBoxToText','on','FontSize',9);
else
    annotation('textbox',dim,'String','No earthquakes recorded','FitBoxToText','on','FontSize',9);
end

