%Input
stn='DAV';                              %abbreviation of the station name, 3 letters. e.g: 'DAV'
station='Davao';                        %full name of the station name. e.g: Davao
country='PHI';                          %abbreviation of country name, 3 letters. e.g: 'IND','PHI'
minmag=5.0;                             %minimum of earthquake magnitude. e.g: 5.0
maxepidis=300;                          %maximum epicentral distance (in km). e.g 150. Insert 'suggest' to get suggestion for maximum epicentral distance as described by Uyeda, 2003

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

dnH=medfilt1(H,'omitnan');
dnD=medfilt1(D,'omitnan');
dnZ=medfilt1(Z,'omitnan');
dnF=medfilt1(F,'omitnan');

sigH=std(dnH,'omitnan');
muH=mean(dnH,'omitnan');
sigD=std(dnD,'omitnan');
muD=mean(dnD,'omitnan');
sigZ=std(dnZ,'omitnan');
muZ=mean(dnZ,'omitnan');
sigF=std(dnF,'omitnan');
muF=mean(dnF,'omitnan');
for i=1:length(dnH)
    if dnH(i)>muH+sigH||dnH(i)<muH-sigH
        dnH(i)=NaN;
    end
    if dnD(i)>muD+sigD||dnD(i)<muD-sigD
        dnD(i)=NaN;
    end
    if dnZ(i)>muZ+sigZ||dnZ(i)<muZ-sigZ
        dnZ(i)=NaN;
    end
    if dnF(i)>muF+sigF||dnF(i)<muF-sigF
        dnF(i)=NaN;
    end
end


%reshape arrays into daily matrices
dayH=reshape(H,86400,[]);
dayD=reshape(D,86400,[]);
dayZ=reshape(Z,86400,[]);
dayF=reshape(F,86400,[]);

%nighttime data only 2300 - 0200 hrs LT 
tStart=0;
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

%band-pass filtering
filtN=2;
filtfL=0.01;
filtfH=0.1;

[b,a]=butter(filtN,[filtfL filtfH]/0.5,'bandpass');

fH=filter(b,a,nightH);
fD=filter(b,a,nightD);
fZ=filter(b,a,nightZ);
fF=filter(b,a,nightF); 

%daily mean
% daymeanH=mean(fH,1,'omitnan');
% daymeanD=mean(fD,1,'omitnan');
% daymeanZ=mean(fZ,1,'omitnan');
% daymeanF=mean(fF,1,'omitnan');

daymeanH=mean(nightH,1,'omitnan');
daymeanD=mean(nightD,1,'omitnan');
daymeanZ=mean(nightZ,1,'omitnan');
daymeanF=mean(nightF,1,'omitnan');

%reshape arrays into whole period data
wholeH=reshape(nightH,[],1);
wholeD=reshape(nightD,[],1);
wholeZ=reshape(nightZ,[],1);
wholeF=reshape(nightF,[],1);

%mean for whole period
wholemeanH=mean(wholeH,1,'omitnan');
wholemeanD=mean(wholeD,1,'omitnan');
wholemeanZ=mean(wholeZ,1,'omitnan');
wholemeanF=mean(wholeF,1,'omitnan');

%standard deviation for whole period
wholestdH=std(wholeH,0,1,'omitnan');
wholestdD=std(wholeD,0,1,'omitnan');
wholestdZ=std(wholeZ,0,1,'omitnan');
wholestdF=std(wholeF,0,1,'omitnan');

%normalisation
normH=(daymeanH-wholemeanH)/wholestdH;
normD=(daymeanD-wholemeanD)/wholestdD;
normZ=(daymeanZ-wholemeanZ)/wholestdZ;
normF=(daymeanF-wholemeanF)/wholestdF;

%H and D exceeding 0.1 only and polarisation ratio
ZH=NaN(1,numdays);
ZD=NaN(1,numdays);

for i=1:numdays
    
    if normH(1,i)>0.1
        ZH(1,i)=normZ(1,i)/normH(1,i);
    end
    if normH(1,i)<=0.1
        ZH(1,i)=0;
    end

    if normD(1,i)>0.1
        ZD(1,i)=normZ(1,i)/normD(1,i);
    end
    if normD(1,i)<=0.1
        ZD(1,i)=0;
    end

end


%mu+/-2sigma
totmeanZH=mean(ZH,'omitnan');
totmeanZD=mean(ZD,'omitnan');
totstdZH=std(ZH,'omitnan');
totstdZD=std(ZD,'omitnan');
mu2psigZH=ones(size(days))*(totmeanZH+2*totstdZH);
mu2msigZH=ones(size(days))*(totmeanZH-2*totstdZH);
mu2psigZD=ones(size(days))*(totmeanZD+2*totstdZD);
mu2msigZD=ones(size(days))*(totmeanZD-2*totstdZD);

%plotting
subplot(3,1,1);
titlestr = sprintf('%s station, %d | Earthquakes: Mw_m_i_n=%.1f, d_m_a_x=%.0f km',station,year,minmag,maxepidis);
plot(days,ZH,'r',days,mu2psigZH,'k--',days,mu2msigZH,'k--');
axis([0 numdays -inf inf]);
set(gca,'XTick',[dataeq(:,2)],'XTickLabel',{eqlbl});
if (min(mu2msigZH)<0 && min(mu2psigZH)>0)
    set(gca,'YTick',[min(mu2msigZH) 0 min(mu2psigZH)],'YTickLabel',{'\mu-2\sigma','0','\mu+2\sigma'});
    else if (min(mu2msigZD)<0 && min(mu2psigZD)<0)
        set(gca,'YTick',[min(mu2msigZH) min(mu2psigZH) 0],'YTickLabel',{'\mu-2\sigma','\mu+2\sigma','0'});
        else set(gca,'YTick',[0 min(mu2msigZD) min(mu2psigZD)],'YTickLabel',{'0','\mu-2\sigma','\mu+2\sigma'});
    end
end   
ylabel('^{Z}/_{H}');
title(titlestr);

subplot(3,1,2);
plot(days,ZD,'r',days,mu2psigZD,'k--',days,mu2msigZD,'k--');
axis([0 numdays -inf inf]);
set(gca,'XTick',[dataeq(:,2)],'XTickLabel',{eqlbl});
if (min(mu2msigZD)<0 && min(mu2psigZD)>0)
    set(gca,'YTick',[min(mu2msigZD) 0 min(mu2psigZD)],'YTickLabel',{'\mu-2\sigma','0','\mu+2\sigma'});
    else if (min(mu2msigZD)<0 && min(mu2psigZD)<0)
        set(gca,'YTick',[min(mu2msigZD) min(mu2psigZD) 0],'YTickLabel',{'\mu-2\sigma','\mu+2\sigma','0'});
        else set(gca,'YTick',[0 min(mu2msigZD) min(mu2psigZD)],'YTickLabel',{'0','\mu-2\sigma','\mu+2\sigma'});
    end
end        
ylabel('^{Z}/_{D}');

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

subplot(3,1,3);
bar(days,Kp,'blue');
hold on
bar(days,Dst,'red');
hold off
xlabel('Days');
ylabel('Dst             \SigmaK_{p}');
legend('\SigmaK_{p}','Dst','location','southeast');
axis([0 numdays -inf inf]);

Details=eqstring;
dim = [0.905 0.725 0.2 0.2];
if numel(Details)~=0
    annotation('textbox',dim,'String',Details,'FitBoxToText','on','FontSize',9);
else
    annotation('textbox',dim,'String','No earthquakes recorded','FitBoxToText','on','FontSize',9);
end


