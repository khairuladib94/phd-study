%Input
stn='DAV';                              %abbreviation of the station name, 3 letters. e.g: 'DAV'
station='Davao';                        %full name of the station name. e.g: Davao
rmstn='MUT';                            %abbreviation of remtote station name
remstation='Muntinlupa';                %full name of remote station for comparison
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

%check remote station
if rmstn=='TGG'
    rmstnnum=1;
end
if rmstn=='MUT'
    rmstnnum=2;
end
if rmstn=='LGZ'
    rmstnnum=3;
end
if rmstn=='CEB'
    rmstnnum=4;
end
if rmstn=='CDO'
    rmstnnum=5;
end
if rmstn=='DAV'
    rmstnnum=6;
end
if rmstn=='GSI'
    rmstnnum=7;
end
if rmstn=='SCN'
    rmstnnum=8;
end
if rmstn=='LWA'
    rmstnnum=9;
end
if rmstn=='PTN'
    rmstnnum=10;
end
if rmstn=='MND';
    rmstnnum=11;
end
if rmstn=='BIK'
    rmstnnum=12;
end
if rmstn=='JYP'
    rmstnnum=13;
end
if rmstn=='PRP'
    rmstnnum=14;
end
if rmstn=='KPG'
    rmstnnum=15;
end
    
latlonstn=[StnIndPhillip(stnnum,2) StnIndPhillip(stnnum,3)];
latlonrmstn=[StnIndPhillip(rmstnnum,2) StnIndPhillip(rmstnnum,3)];
disstn=deg2km(distance('gc',latlonstn,latlonrmstn));

%check leap year
if rem(year,4)==0
    days=2:1:366;
    numdays=366-1;
else
    days=2:1:365;
    numdays=365-1;
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
rmdataeq=NaN(100,4);
dataeqnum=1;
rmdataeqnum=1;
for numeq=1:toteq
    if eqvar(numeq,1)==year
        if  eqvar(numeq,10)>=minmag
            latloneq=[eqvar(numeq,7) eqvar(numeq,8)];
            epidis=deg2km(distance('gc',latlonstn,latloneq));
            rmepidis=deg2km(distance('gc',latlonrmstn,latloneq));
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
    
            if rmepidis<=maxepidis
                rmdataeq(rmdataeqnum,1)=dateeq;
                rmdataeq(rmdataeqnum,2)=daynum;
                rmdataeq(rmdataeqnum,3)=rmepidis;
                rmdataeq(rmdataeqnum,4)=mag;
                rmdataeqnum=rmdataeqnum+1;
            end
            
        end
    end
end

dataeq(any(isnan(dataeq),2),:)=[];
if size(dataeq,1)>1
    dataeq=flip(dataeq);
end

rmdataeq(any(isnan(rmdataeq),2),:)=[];
if size(rmdataeq,1)>1
    rmdataeq=flip(rmdataeq);
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

rmeqstring=string(zeros(size(rmdataeq,1),1));
eqdigit=string(1:50);
rmeql=zeros(size(rmdataeq,1),1);
rmeql=(string(eqdigit(1:(size(rmdataeq,1)))))';

eqlbl=eql;
rmeqlbl=rmeql;

for counth=1:size(dataeq,1)
    if counth~=1
     if (dataeq(counth,2)-dataeq(counth-1,2))<=2
        eqlbl(counth)=sprintf('%s%s',eqlbl(counth-1),eqlbl(counth));
        eqlbl(counth-1)=' ';
     end
    end
end

for countp=1:size(rmdataeq,1)
    if countp~=1
     if (rmdataeq(countp,2)-rmdataeq(countp-1,2))<=2
        rmeqlbl(countp)=sprintf('%s%s',eqlbl(countp-1),rmeqlbl(countp));
        eqlbl(countp-1)=' ';
     end
    end
end

for countf=1:size(dataeq,1)
   date=datetime(datevec(dataeq(countf,1)),'Format','dd/MM');
   eqstring(countf,1)=sprintf('%s, M=%.1f, %.0f km',date,dataeq(countf,4),dataeq(countf,3));
end

for counto=1:size(rmdataeq,1)
   date=datetime(datevec(rmdataeq(counto,1)),'Format','dd/MM');
   rmeqstring(counto,1)=sprintf('%s, M=%.1f, %.0f km',date,rmdataeq(counto,4),rmdataeq(counto,3));
end

if numel(eql)>0
    eqstring=strcat('(',eql,')',eqstring);
end

if numel(rmeql)>0
    rmeqstring=strcat('(',rmeql,')',rmeqstring);
end

%----------------------------------------------------------------------------------------------------------------------------------------------
%Geomagnetic field

%remove noise/error

dnH1=medfilt1(H1,'omitnan');
dnD1=medfilt1(D1,'omitnan');
dnZ1=medfilt1(Z1,'omitnan');

dnH2=medfilt1(H2,'omitnan');
dnD2=medfilt1(D2,'omitnan');
dnZ2=medfilt1(Z2,'omitnan');

for k=1:3
    
    sigH1=std(dnH1,'omitnan');
    muH1=mean(dnH1,'omitnan');
    sigD1=std(dnD1,'omitnan');
    muD1=mean(dnD1,'omitnan');
    sigZ1=std(dnZ1,'omitnan');
    muZ1=mean(dnZ1,'omitnan');
    
    sigH2=std(dnH2,'omitnan');
    muH2=mean(dnH2,'omitnan');
    sigD2=std(dnD2,'omitnan');
    muD2=mean(dnD2,'omitnan');
    sigZ2=std(dnZ2,'omitnan');
    muZ2=mean(dnZ2,'omitnan');
    
    for i=1:length(dnH1)
        if dnH1(i)>muH1+4*sigH1||dnH1(i)<muH1-3*sigH1
            dnH1(i)=NaN;
        end
        if dnD1(i)>muD1+4*sigD1||dnD1(i)<muD1-3*sigD1
            dnD1(i)=NaN;
        end
        if dnZ1(i)>muZ1+4*sigZ1||dnZ1(i)<muZ1-3*sigZ1
            dnZ1(i)=NaN;
        end
        if dnH2(i)>muH2+4*sigH2||dnH2(i)<muH2-3*sigH2
            dnH2(i)=NaN;
        end
        if dnD2(i)>muD2+4*sigD2||dnD2(i)<muD2-3*sigD2
            dnD2(i)=NaN;
        end
        if dnZ2(i)>muZ2+4*sigZ2||dnZ2(i)<muZ2-3*sigZ2
            dnZ2(i)=NaN;
        end
    end
    
end

%close station diurnal variation
H_1=dnH1(57601:length(dnH1)-28800);
D_1=dnD1(57601:length(dnD1)-28800);
Z_1=dnZ1(57601:length(dnZ1)-28800);

H_1=reshape(H_1,[],numdays); 
D_1=reshape(D_1,[],numdays);
Z_1=reshape(Z_1,[],numdays);

delH1=NaN(numdays);
delD1=NaN(numdays);
delZ1=NaN(numdays);

H1min=min(H_1);
H1max=max(H_1);
delH1=H1max-H1min;

D1min=min(D_1);
D1max=max(D_1);
delD1=D1max-D1min;

Z1min=min(Z_1);
Z1max=max(Z_1);
delZ1=Z1max-Z1min;

%remote station diurnal variation
H_2=dnH2(57601:length(dnH2)-28800);
D_2=dnD2(57601:length(dnD2)-28800);
Z_2=dnZ2(57601:length(dnZ2)-28800);

H_2=reshape(H_2,[],numdays);
D_2=reshape(D_2,[],numdays);
Z_2=reshape(Z_2,[],numdays);

delH2=NaN(numdays);
delD2=NaN(numdays);
delZ2=NaN(numdays);

H2min=min(H_2);
H2max=max(H_2);
delH2=H2max-H2min;

D2min=min(D_2);
D2max=max(D_2);
delD2=D2max-D2min;

Z2min=min(Z_2);
Z2max=max(Z_2);
delZ2=Z2max-Z2min;

%ratio
H1H2=delH1./delH2;
D1D2=delD1./delD2;
Z1Z2=delZ1./delZ2;

mmdelH1=movmean(delH1,10,'includenan');
mmdelD1=movmean(delD1,10,'includenan');
mmdelZ1=movmean(delZ1,10,'includenan');

mmdelH2=movmean(delH2,10,'includenan');
mmdelD2=movmean(delD2,10,'includenan');
mmdelZ2=movmean(delZ2,10,'includenan');

mmH1H2=movmean(H1H2,10,'includenan');
mmD1D2=movmean(D1D2,10,'includenan');
mmZ1Z2=movmean(Z1Z2,10,'includenan');

%mu+/-3sigma
muH1H2=mean(H1H2,'omitnan');
muD1D2=mean(D1D2,'omitnan');
muZ1Z2=mean(Z1Z2,'omitnan');
sigH1H2=std(H1H2,'omitnan');
sigD1D2=std(D1D2,'omitnan');
sigZ1Z2=std(Z1Z2,'omitnan');

mup3sigH1H2=ones(size(days))*(muH1H2+3*sigH1H2);
mum3sigH1H2=ones(size(days))*(muH1H2-3*sigH1H2);
mup3sigD1D2=ones(size(days))*(muD1D2+3*sigD1D2);
mum3sigD1D2=ones(size(days))*(muD1D2-3*sigD1D2);
mup3sigZ1Z2=ones(size(days))*(muZ1Z2+3*sigZ1Z2);
mum3sigZ1Z2=ones(size(days))*(muZ1Z2-3*sigZ1Z2);

%----------------------------------------------------------------------------------------------------------------------------------------------
%Plotting

figure(1)

subplot(4,1,1);
titlestr = sprintf('%s station, %d. Remote station: %s (%.0f km away) | Earthquakes: Mw_m_i_n=%.1f, d_m_a_x=%.0f km',station,year,remstation,disstn,minmag,maxepidis);
plot(days,H1H2,days,mup3sigH1H2,'k--',days,mum3sigH1H2,'k--');
hold on
plot(days,mmH1H2,'r');
hold off
axis([0 numdays -inf inf]);
set(gca,'XTick',[dataeq(:,2)],'XTickLabel',{eqlbl});
ylabel('\DeltaH_{near}/\DeltaH_{remote}');
if (min(mum3sigH1H2)<0 && min(mup3sigH1H2)>0)
    set(gca,'YTick',[min(mum3sigH1H2) 0 min(mup3sigH1H2)],'YTickLabel',{'\mu-3\sigma','0','\mu+3\sigma'});
    else if (min(mum3sigH1H2)<0 && min(mup3sigH1H2)<0)
        set(gca,'YTick',[min(mum3sigH1H2) min(mup3sigH1H2) 0],'YTickLabel',{'\mu-3\sigma','\mu+3\sigma','0'});
        else set(gca,'YTick',[0 min(mum3sigH1H2) min(mup3sigH1H2)],'YTickLabel',{'0','\mu-3\sigma','\mu+3\sigma'});
    end
end   
title(titlestr);

subplot(4,1,2);
plot(days,D1D2,days,mup3sigD1D2,'k--',days,mum3sigD1D2,'k--');
hold on
plot(days,mmD1D2,'r');
hold off
axis([0 numdays -inf inf]);
set(gca,'XTick',[dataeq(:,2)],'XTickLabel',{eqlbl});
ylabel('\DeltaD_{near}/\DeltaD_{remote}');
if (min(mum3sigD1D2)<0 && min(mup3sigD1D2)>0)
    set(gca,'YTick',[min(mum3sigD1D2) 0 min(mup3sigD1D2)],'YTickLabel',{'\mu-3\sigma','0','\mu+3\sigma'});
    else if (min(mum3sigD1D2)<0 && min(mup3sigD1D2)<0)
        set(gca,'YTick',[min(mum3sigD1D2) min(mup3sigD1D2) 0],'YTickLabel',{'\mu-3\sigma','\mu+3\sigma','0'});
        else set(gca,'YTick',[0 min(mum3sigD1D2) min(mup3sigD1D2)],'YTickLabel',{'0','\mu-3\sigma','\mu+3\sigma'});
    end
end   

subplot(4,1,3);
plot(days,Z1Z2,days,mup3sigZ1Z2,'k--',days,mum3sigZ1Z2,'k--');
hold on
plot(days,mmZ1Z2,'r');
hold off
axis([0 numdays -inf inf]);
set(gca,'XTick',[dataeq(:,2)],'XTickLabel',{eqlbl});
ylabel('\DeltaZ_{near}/\DeltaZ_{remote}');
if (min(mum3sigZ1Z2)<0 && min(mup3sigZ1Z2)>0)
    set(gca,'YTick',[min(mum3sigZ1Z2) 0 min(mup3sigZ1Z2)],'YTickLabel',{'\mu-3\sigma','0','\mu+3\sigma'});
    else if (min(mum3sigZ1Z2)<0 && min(mup3sigZ1Z2)<0)
        set(gca,'YTick',[min(mum3sigZ1Z2) min(mup3sigZ1Z2) 0],'YTickLabel',{'\mu-3\sigma','\mu+3\sigma','0'});
        else set(gca,'YTick',[0 min(mum3sigZ1Z2) min(mup3sigZ1Z2)],'YTickLabel',{'0','\mu-3\sigma','\mu+3\sigma'});
    end
end   

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
bar(1:365,Kp,'blue');
hold on
bar(1:365,Dst,'red');
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

rmDetails=rmeqstring;
dim = [0.905 0.1 0.2 0.2];
if numel(rmDetails)~=0
    annotation('textbox',dim,'String',rmDetails,'FitBoxToText','on','FontSize',9);
else
    annotation('textbox',dim,'String','No earthquakes recorded','FitBoxToText','on','FontSize',9);
end


figure(2)

subplot(4,1,1);
titlestr = sprintf('%s station, %d. Remote station: %s (%.0f km away) | Earthquakes: Mw_m_i_n=%.1f, d_m_a_x=%.0f km',station,year,remstation,disstn,minmag,maxepidis);
plot(days,delH1);
hold on
plot(days,delH2,'r');
hold off
axis([0 numdays -inf inf]);
set(gca,'XTick',[dataeq(:,2)],'XTickLabel',{eqlbl});
ylabel('\DeltaH');
legend(sprintf('%s',stn),sprintf('%s',rmstn),'location','northeast');
title(titlestr);

subplot(4,1,2);
plot(days,delD1);
hold on
plot(days,delD2,'r');
hold off
axis([0 numdays -inf inf]);
set(gca,'XTick',[dataeq(:,2)],'XTickLabel',{eqlbl});
ylabel('\DeltaD');
legend(sprintf('%s',stn),sprintf('%s',rmstn),'location','northeast');

subplot(4,1,3);
plot(days,delZ1);
hold on
plot(days,delZ2,'r');
hold off
axis([0 numdays -inf inf]);
set(gca,'XTick',[dataeq(:,2)],'XTickLabel',{eqlbl});
ylabel('\DeltaZ');
legend(sprintf('%s',stn),sprintf('%s',rmstn),'location','northeast');

Kp=zeros(numdays+1,1);
Dst=zeros(numdays+1,1);
counte=1;
for countd=1:length(Ind20112017)    
    if Ind20112017(countd,1)==year
         Kp(counte,1)=Ind20112017(countd,5);
         Dst(counte,1)=Ind20112017(countd,6);
         counte=counte+1;
    end  
end

subplot(4,1,4);
bar(1:365,Kp,'blue');
hold on
bar(1:365,Dst,'red');
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

rmDetails=rmeqstring;
dim = [0.905 0.05 0.2 0.2];
if numel(rmDetails)~=0
    annotation('textbox',dim,'String',rmDetails,'FitBoxToText','on','FontSize',9);
else
    annotation('textbox',dim,'String','No earthquakes recorded','FitBoxToText','on','FontSize',9);
end

