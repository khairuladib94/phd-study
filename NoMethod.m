%Input
stn='CDO';                              %abbreviation of the station name, 3 letters. e.g: 'DAV'
station='Cagayan De Oro';                        %full name of the station name. e.g: Davao
country='PHI';                          %abbreviation of country name, 3 letters. e.g: 'IND','PHI'
minmag=7.0;                             %minimum of earthquake magnitude. e.g: 5.0
maxepidis=1000;                          %maximum epicentral distance (in km). e.g 150. Insert 'suggest' to get suggestion for maximum epicentral distance as described by Uyeda, 2003

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
    secs=0:1:86400*366-1;
    numsecs=86400*366;
    days=1:1:366;
    numdays=366;
else
    secs=0:1:86400*365-1;
    numsecs=86400*365;
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

dataeq=NaN(100,4);
dataeqnum=1;
for numeq=1:length(eqvar)
    if eqvar(numeq,1)==year
        if  eqvar(numeq,10)>=minmag
            latloneq=[eqvar(numeq,7) eqvar(numeq,8)];
            epidis=deg2km(distance('gc',latlonstn,latloneq));
            dateeq=datenum([eqvar(numeq,1) eqvar(numeq,2) eqvar(numeq,3)]);
            mag=eqvar(numeq,10);
            eqday=day(datetime([eqvar(numeq,1) eqvar(numeq,2) eqvar(numeq,3)]),'dayofyear');
            eqhr=datenum(eqvar(numeq,4));
            eqmin=datenum(eqvar(numeq,5));
            eqsec=datenum(eqvar(numeq,6));
            eqsecofyear=((eqday-1)*86400)+(eqhr*3600)+(eqmin*60)+eqsec;
            
            if epidis<=maxepidis
                dataeq(dataeqnum,1)=dateeq;
                dataeq(dataeqnum,2)=eqsecofyear;
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

%remove noise/error
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

%band-pass filtering
filtN=3;
filtfL=0.01;
filtfH=0.1;

[b,a]=butter(filtN,[filtfL filtfH]/0.5,'bandpass');

% dnH=filter(b,a,dnH);
% dnD=filter(b,a,dnD);
% dnZ=filter(b,a,dnZ);
% dnF=filter(b,a,dnF); 

%plotting
subplot(5,1,1);
titlestr = sprintf('%s station, %d | Earthquakes: Mw_m_i_n=%.1f, d_m_a_x=%.0f km',station,year,minmag,maxepidis);
plot(secs,dnH);
axis([0 numsecs -inf inf]);
set(gca,'XTick',[dataeq(:,2)],'XTickLabel',{eqlbl});
ylabel('H');
title(titlestr);

subplot(5,1,2);
plot(secs,dnD);
axis([0 numsecs -inf inf]);
set(gca,'XTick',[dataeq(:,2)],'XTickLabel',{eqlbl});
ylabel('D');

subplot(5,1,3);
plot(secs,dnZ);
axis([0 numsecs -inf inf]);
set(gca,'XTick',[dataeq(:,2)],'XTickLabel',{eqlbl});
ylabel('Z');

subplot(5,1,4);
plot(secs,dnF);
axis([0 numsecs -inf inf]);
set(gca,'XTick',[dataeq(:,2)],'XTickLabel',{eqlbl});
ylabel('F');

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

subplot(5,1,5);
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

