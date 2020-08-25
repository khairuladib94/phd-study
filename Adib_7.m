%Place and time inputs-----------------------------------------------------

stn='LWA';                       %Abbreviation of station name
date_start=[2014,01,01];           %Insert custom start and end dates
date_end=[2014,05,31];           %Period spanning through 3 consecutive years is the maximum  

%Customization-------------------------------------------------------------

mag_min=5.0;                     %Minimum magnitude of earthquakes to be considered
dis_max=300;                     %Maximum epicentral distance from the station

nfft=2048;
size_win=1800;
noverlap=0;

f_A1=0.01;
f_A2=0.1;

f_B1=0.01;
f_B2=0.02;

f_C1=0.09;
f_C2=0.1;

%Time period---------------------------------------------------------------

datenum_start=datenum(date_start);
datenum_end=datenum(date_end);
if date_start(1)==date_end(1)
    year=string(date_start(1));
else
    year=sprintf('%s-%s',string(date_start(1)),string(date_end(1)));
end

%Loading files-------------------------------------------------------------

load VARIABLES_WORLD

j=date_start(1);
for i=1:3
    year_vec(1,i)=j;
    if j==date_end(1)
        break
    end
    j=j+1;
end

if numel(year_vec)==1
    matname=strcat(stn,string(year_vec(1)),'S');
    load(matname);
end
    
if numel(year_vec)==2
    matname(1)=strcat(stn,string(year_vec(1)),'S');
    matname(2)=strcat(stn,string(year_vec(2)),'S');
    A=load(matname(1));
    B=load(matname(2));
    H=vertcat(A.H,B.H);
    D=vertcat(A.D,B.D);
    Z=vertcat(A.Z,B.Z);
    UT1m=horzcat(A.UT1m,B.UT1m);
end

if numel(year_vec)==3
    matname(1)=strcat(stn,string(year_vec(1)),'S');
    matname(2)=strcat(stn,string(year_vec(2)),'S');
    matname(3)=strcat(stn,string(year_vec(3)),'S');
    A=load(matname(1));
    B=load(matname(2));
    C=load(matname(3));
    H=vertcat(A.H,B.H,C.H);
    D=vertcat(A.D,B.D,C.D);
    Z=vertcat(A.Z,B.Z,C.Z);
    UT1m=horzcat(A.UT1m,B.UT1m,C.UT1m);
end

%--------------------------------------------------------------------------

days_num=datenum_end-datenum_start+1;
days_vec=datetime(datevec(datenum_start:1:datenum_end),'Format','dd/MM/yyyy');

%Setting station number

stn_vec={'TGG';'MUT';'LGZ';'CEB';'CDO';'DAV';'GSI';'SCN';'LWA';'PTN';'MND';
    'BIK';'JYP';'PRP';'KPG';'KTB';'DAW';'LKW';'SBH';'PER';'BTG';'KTN';'TIK';
    'CHD';'CST';'ZYK';'ZGN';'MGD';'YAK';'PTK';'ASB';'TNO';'ONW';'OIS';'KUJ';
    'AMA';'HLN';'EWA';'YAP';'BCL';'HVD';'TIR';'CMB';'CKT';'TWV';'ROC';'LMT';
    'CGR';'CMD';'CAN';'MLB';'HOB';'MCQ';'DVS';'WAD';'GLY';'JRS';'TPT';'TMA';
    'ANC';'HUA';'ICA';'EUS';'SMA';'LAQ';'FYM';'ASW';'KRT';'AAB';'ILR';'ABU';
    'LAG';'ABJ';'NAB';'DES';'LSK';'MPT';'DRB';'HER'};

station_vec={'Tuguegarao';'Muntinlupa';'Legazpi';'Cebu';'Cagayan De Oro';'Davao';
    'Gunung Sitoli';'Sicincin';'Liwa';'Pontianak';'Manado';'Biak';'Jayapura';'Pare Pare';
    'Kupang';'Kototabang';'Darwin';'Langkawi';'Sabah';'Perak';'Banting';'Kotelnyy';
    'Tixie';'Chokurdakh';'Cape Schmidt';'Zyryanka';'Zhigansk';'Magadan';'Yakutsk';
    'Paratunka';'Ashibetsu';'Tohno';'Onagawa';'Oiso';'Kuju';'Amami-Oh-shima';
    'Hualien';'Ewa Beach';'Yap Island';'Bac Lieu';'Khovd';'Tirunelveli';'Colombo';
    'Cooktown';'Townsville';'Rockhampton';'Learmonth';'Culgoora';'Camden';'Canberra';
    'Crib Point';'Hobart';'Macquarie Island';'Davis';'Wadena';'Glyndon';'Jerusalem';
    'Tarapoto';'Tingo Maria';'Ancon';'Huancayo';'Ica';'Eusebio';'Santa Maria';
    'Laquila';'Fayum';'Aswan';'Khartoum';'Adis Ababa';'Ilorin';'Abuja';'Lagos';
    'Abidjan';'Nairobi';'Dal Es Salaam';'Lusaka';'Maputo';'Durban';'Hermanus'};

for i=1:length(stn_MAGDAS)
    if strcmp(stn,stn_vec(i))
        stn_num=i;
        station=string(station_vec(i));
        break;
    end
end
    
%Getting station coordinate
stn_latlon=[stn_MAGDAS(stn_num,2:3)];

%Earthquake----------------------------------------------------------------

%Building selected earthquakes table

j=1;
for i=1:length(EQ_WORLD)
    if datenum(EQ_WORLD(i,1:3))>=datenum_start && datenum(EQ_WORLD(i,1:3))<=datenum_end
        EQ_table(j,:)=EQ_WORLD(i,:);
        j=j+1;
    end
    if datenum(EQ_WORLD(i,1:3))>datenum_end
        break
    end
end

i=1;
for j=1:length(EQ_table)
    if  EQ_table(j,10)>=mag_min
        EQ_latlon=[EQ_table(j,7:8)];
        EQ_dis=deg2km(distance('gc',stn_latlon,EQ_latlon));
        EQ_date=datenum(EQ_table(j,1:3));
        EQ_mag=EQ_table(j,10);
        EQ_DOY=day(datetime([EQ_table(j,1:3)]),'dayofyear');
        %EQ_Ks=((10^(0.75*EQ_mag))/(10*EQ_dis))*((1+EQ_dis*10^(-EQ_mag/2))^(-2.33));
        EQ_Ks=(10^(0.75*EQ_mag))/(EQ_dis+100);
        EQ_depth=EQ_table(j,9);
        EQ_f=1000/(pi*((EQ_depth*1000)^2)*1.2566e-06);
        
        
        if EQ_dis<=dis_max
            EQ_sel(i,1)=EQ_date;
            EQ_sel(i,2)=EQ_DOY;
            EQ_sel(i,3)=EQ_dis;
            EQ_sel(i,4)=EQ_mag;
            EQ_sel(i,5)=EQ_Ks;
            EQ_sel(i,6)=EQ_depth;
            EQ_sel(i,7)=EQ_f;
            EQ_sel(i,8:9)=EQ_latlon;
            i=i+1;
        end
    end
    
end

%Geomagnetic data----------------------------------------------------------

data_start=round(min(UT1m));
data_end=floor(max(UT1m));

day_start=datenum_start-data_start+1;
day_end=datenum_end-data_start+1;
datenum_vec=datenum_start:datenum_end;

H_period=H(86400*day_start-86399:86400*day_end);
D_period=D(86400*day_start-86399:86400*day_end);
Z_period=Z(86400*day_start-86399:86400*day_end);

%Removing noise/outlier
H_dn=medfilt1(H_period,'omitnan');
D_dn=medfilt1(D_period,'omitnan');
Z_dn=medfilt1(Z_period,'omitnan');

for i=1:5
    
    H_sig=std(H_dn,'omitnan');
    H_mu=mean(H_dn,'omitnan');
    D_sig=std(D_dn,'omitnan');
    D_mu=mean(D_dn,'omitnan');
    Z_sig=std(Z_dn,'omitnan');
    Z_mu=mean(Z_dn,'omitnan');
    
    for j=1:length(H_dn)
        
        if H_dn(j)>H_mu+5*H_sig||H_dn(j)<H_mu-5*H_sig
            H_dn(j)=NaN;
        end
        if D_dn(j)>D_mu+5*D_sig||D_dn(j)<D_mu-5*D_sig
            D_dn(j)=NaN;
        end
        if Z_dn(j)>Z_mu+5*Z_sig||Z_dn(j)<Z_mu-5*Z_sig
            Z_dn(j)=NaN;
        end
        
    end
    
end

%Nighttime data 2200 - 0200 hrs LT
t_zone=timezone(stn_latlon(2),'degrees');
t_start=22+t_zone;      
if t_start>24
    t_start=t_start-24;
end
t_start=t_start*3600;      

for i=1:days_num
    
    H_night(:,i)=H_dn(t_start:t_start+14400-1);
    D_night(:,i)=D_dn(t_start:t_start+14400-1);
    Z_night(:,i)=Z_dn(t_start:t_start+14400-1);
    t_start=t_start+86400;
    
end

%---------------------------------------------------------

bin_1=round(((nfft/2+1)/0.5)*f_A1);
bin_2=round(((nfft/2+1)/0.5)*f_A2);

f_s=1;

[H_welch,f]=pwelch(H_night,size_win,noverlap,nfft,f_s);
H_welch=10*log10(real(H_welch));
H_bin=H_welch(bin_1:bin_2,:);
H_daymu=mean(H_bin,1,'omitnan');

[D_welch,f]=pwelch(D_night,size_win,noverlap,nfft,f_s);
D_welch=10*log10(real(D_welch));
D_bin=D_welch(bin_1:bin_2,:);
D_daymu=mean(D_bin,1,'omitnan');

[Z_welch,f]=pwelch(Z_night,size_win,noverlap,nfft,f_s);
Z_welch=10*log10(real(Z_welch));
Z_bin=Z_welch(bin_1:bin_2,:);
Z_daymu=mean(Z_bin,1,'omitnan');

%---------------------------------

%Normalisation

j=0;
k=1;
for i=1:days_num
    m=datevec(datenum_vec(i));
    if i<days_num
        n=datevec(datenum_vec(i+1));
    end
    if i<days_num && m(2)==n(2)
        j=j+1;
    end
    if m(2)~=n(2) || i==days_num
        H_monmu(1,k)=mean(H_daymu(i-j:i),'omitnan');
        D_monmu(1,k)=mean(D_daymu(i-j:i),'omitnan');
        Z_monmu(1,k)=mean(Z_daymu(i-j:i),'omitnan');
        
        H_monsig(1,k)=std(H_daymu(i-j:i),'omitnan');
        D_monsig(1,k)=std(D_daymu(i-j:i),'omitnan');
        Z_monsig(1,k)=std(Z_daymu(i-j:i),'omitnan');
        
        H1_norm(1,i-j:i)=(H_daymu(1,i-j:i)-H_monmu(1,k))/H_monsig(1,k);
        D1_norm(1,i-j:i)=(D_daymu(1,i-j:i)-D_monmu(1,k))/D_monsig(1,k);
        Z1_norm(1,i-j:i)=(Z_daymu(1,i-j:i)-Z_monmu(1,k))/Z_monsig(1,k);
        
        k=k+1;
        j=0;
    end
    
end

%----------------------------------


ZG1=Z1_norm./(sqrt(H1_norm.^2+D1_norm.^2));

%mu+/-2sigma
ZG1_mu=mean(ZG1,'omitnan');
ZG1_sig=std(ZG1,'omitnan');
ZG1_mup2sig=ones(1,days_num)*(ZG1_mu+2*ZG1_sig);
ZG1_mum2sig=ones(1,days_num)*(ZG1_mu-2*ZG1_sig);
ZG1_mm=movmean(ZG1,5,'omitnan');

%------------------------------------

bin_1=round(((nfft/2+1)/0.5)*f_B1);
bin_2=round(((nfft/2+1)/0.5)*f_B2);

[H_welch,f]=pwelch(H_night,size_win,noverlap,nfft,f_s);
H_welch=10*log10(real(H_welch));
H_bin=H_welch(bin_1:bin_2,:);
H_daymu=mean(H_bin,1,'omitnan');

[D_welch,f]=pwelch(D_night,size_win,noverlap,nfft,f_s);
D_welch=10*log10(real(D_welch));
D_bin=D_welch(bin_1:bin_2,:);
D_daymu=mean(D_bin,1,'omitnan');

[Z_welch,f]=pwelch(Z_night,size_win,noverlap,nfft,f_s);
Z_welch=10*log10(real(Z_welch));
Z_bin=Z_welch(bin_1:bin_2,:);
Z_daymu=mean(Z_bin,1,'omitnan');

%---------------------------------

%Normalisation

j=0;
k=1;
for i=1:days_num
    m=datevec(datenum_vec(i));
    if i<days_num
        n=datevec(datenum_vec(i+1));
    end
    if i<days_num && m(2)==n(2)
        j=j+1;
    end
    if m(2)~=n(2) || i==days_num
        H_monmu(1,k)=mean(H_daymu(i-j:i),'omitnan');
        D_monmu(1,k)=mean(D_daymu(i-j:i),'omitnan');
        Z_monmu(1,k)=mean(Z_daymu(i-j:i),'omitnan');
        
        H_monsig(1,k)=std(H_daymu(i-j:i),'omitnan');
        D_monsig(1,k)=std(D_daymu(i-j:i),'omitnan');
        Z_monsig(1,k)=std(Z_daymu(i-j:i),'omitnan');
        
        H2_norm(1,i-j:i)=(H_daymu(1,i-j:i)-H_monmu(1,k))/H_monsig(1,k);
        D2_norm(1,i-j:i)=(D_daymu(1,i-j:i)-D_monmu(1,k))/D_monsig(1,k);
        Z2_norm(1,i-j:i)=(Z_daymu(1,i-j:i)-Z_monmu(1,k))/Z_monsig(1,k);
        
        k=k+1;
        j=0;
    end
    
end

ZG2=Z2_norm./(sqrt(H2_norm.^2+D2_norm.^2));

%mu+/-2sigma
ZG2_mu=mean(ZG2,'omitnan');
ZG2_sig=std(ZG2,'omitnan');
ZG2_mup2sig=ones(1,days_num)*(ZG2_mu+2*ZG2_sig);
ZG2_mum2sig=ones(1,days_num)*(ZG2_mu-2*ZG2_sig);
ZG2_mm=movmean(ZG2,5,'omitnan');


%----------------------------------------

bin_1=round(((nfft/2+1)/0.5)*f_C1);
bin_2=round(((nfft/2+1)/0.5)*f_C2);

[H_welch,f]=pwelch(H_night,size_win,noverlap,nfft,f_s);
H_welch=10*log10(real(H_welch));
H_bin=H_welch(bin_1:bin_2,:);
H_daymu=mean(H_bin,1,'omitnan');

[D_welch,f]=pwelch(D_night,size_win,noverlap,nfft,f_s);
D_welch=10*log10(real(D_welch));
D_bin=D_welch(bin_1:bin_2,:);
D_daymu=mean(D_bin,1,'omitnan');

[Z_welch,f]=pwelch(Z_night,size_win,noverlap,nfft,f_s);
Z_welch=10*log10(real(Z_welch));
Z_bin=Z_welch(bin_1:bin_2,:);
Z_daymu=mean(Z_bin,1,'omitnan');

%---------------------------------

%Normalisation

j=0;
k=1;
for i=1:days_num
    m=datevec(datenum_vec(i));
    if i<days_num
        n=datevec(datenum_vec(i+1));
    end
    if i<days_num && m(2)==n(2)
        j=j+1;
    end
    if m(2)~=n(2) || i==days_num
        H_monmu(1,k)=mean(H_daymu(i-j:i),'omitnan');
        D_monmu(1,k)=mean(D_daymu(i-j:i),'omitnan');
        Z_monmu(1,k)=mean(Z_daymu(i-j:i),'omitnan');
        
        H_monsig(1,k)=std(H_daymu(i-j:i),'omitnan');
        D_monsig(1,k)=std(D_daymu(i-j:i),'omitnan');
        Z_monsig(1,k)=std(Z_daymu(i-j:i),'omitnan');
        
        H3_norm(1,i-j:i)=(H_daymu(1,i-j:i)-H_monmu(1,k))/H_monsig(1,k);
        D3_norm(1,i-j:i)=(D_daymu(1,i-j:i)-D_monmu(1,k))/D_monsig(1,k);
        Z3_norm(1,i-j:i)=(Z_daymu(1,i-j:i)-Z_monmu(1,k))/Z_monsig(1,k);
        
        k=k+1;
        j=0;
    end
    
end

ZG3=Z3_norm./(sqrt(H3_norm.^2+D3_norm.^2));

%mu+/-2sigma
ZG3_mu=mean(ZG3,'omitnan');
ZG3_sig=std(ZG3,'omitnan');
ZG3_mup2sig=ones(1,days_num)*(ZG3_mu+2*ZG3_sig);
ZG3_mum2sig=ones(1,days_num)*(ZG3_mu-2*ZG3_sig);
ZG3_mm=movmean(ZG3,5,'omitnan');


%----------------------------------------------

%Plotting seismic index chart

f1=figure(1);

subplot(5,1,1)
for i=1:size(EQ_sel,1)
    j=datetime(datevec(EQ_sel(i,1)),'Format','dd/MM/yyyy');
    if EQ_sel(i,6)<=50
        plot(j,EQ_sel(i,5),'ro');
    elseif (EQ_sel(i,6)>50) && (EQ_sel(i,6)<=200)
        plot(j,EQ_sel(i,5),'b^');
    else EQ_sel(i,6)>200
         plot(j,EQ_sel(i,5),'gx');
    end
    
    if i==1
       hold on
    end
end
xlim([min(days_vec) max(days_vec)])
set(gca, 'XTickLabel',[]);
ylim([0 max(EQ_sel(:,5))+10]);
ylabel('K_{s}');
h=zeros(3,1);
h(1) = plot(days_vec,NaN(1,days_num),'ro');
h(2) = plot(days_vec,NaN(1,days_num),'b^');
h(3) = plot(days_vec,NaN(1,days_num),'gx');
lgd=legend(h,'<50km','50-200km','>200km','location','northwest');
title(lgd,'Depth');
title(sprintf('^{Z}/_{G} Polarisation ratio | %s station, %s - %s',station,min(days_vec),max(days_vec)));
hold off


%Plotting polarisation ratio

subplot(5,1,2);
plot(days_vec,ZG1,days_vec,ZG1_mum2sig,'r--',days_vec,ZG1_mup2sig,'r--');
ylabel('^{Z}/_{G}');
legend(sprintf('\\Deltaf=%.3f-%.3f Hz',f_A1,f_A2),'location','northwest');
xlim([min(days_vec) max(days_vec)])
set(gca, 'XTickLabel',[]);

subplot(5,1,3);
plot(days_vec,ZG2,days_vec,ZG2_mum2sig,'r--',days_vec,ZG2_mup2sig,'r--');
ylabel('^{Z}/_{G}');
legend(sprintf('\\Deltaf=%.3f-%.3f Hz',f_B1,f_B2),'location','northwest');
xlim([min(days_vec) max(days_vec)])
set(gca, 'XTickLabel',[]);

subplot(5,1,4);
plot(days_vec,ZG3,days_vec,ZG3_mum2sig,'r--',days_vec,ZG3_mup2sig,'r--');
ylabel('^{Z}/_{G}');
legend(sprintf('\\Deltaf=%.3f-%.3f Hz',f_C1,f_C2),'location','northwest');
xlim([min(days_vec) max(days_vec)])
set(gca, 'XTickLabel',[]);

%Getting Kp and Dst indices
j=1;
for i=1:length(index_geomag)
    m=datenum(index_geomag(i,1),1,index_geomag(i,2));
    if m>=datenum_start && m<=datenum_end
        Kp(j,1)=index_geomag(i,4);
        Dst(j,1)=index_geomag(i,5);
        j=j+1;
        
        if j==days_num+1
            break;
        end
    end
end

%Plotting Kp and Dst indices

subplot(5,1,5)
bar(days_vec,Kp,'b');
hold on
bar(days_vec,Dst,'r');
plot(days_vec,ones(days_num,1)*(32),'b--');
plot(days_vec,ones(days_num,1)*(-50),'r--');
hold off
xlabel('Date');
ylabel('Dst     \SigmaK_{p}');
legend('\SigmaK_{p}','Dst','location','northwest');
xlim([min(days_vec) max(days_vec)])
ylim([-60 40]);


% plot(days_vec,Z1_norm,days_vec,1./(sqrt(H1_norm.^2+D1_norm.^2)),'r');
% xlim([min(days_vec) max(days_vec)])
% legend('Z','^{1}/_{G}','location','northwest');
% title(sprintf('Normalized Z and ^{1}/_{G} plots | %s station, %s - %s',station,min(days_vec),max(days_vec)));
% xlabel('Date');

set(gcf, 'Position', get(0, 'Screensize'));

%--------------------------------------------------------------------------


%Create map

f2=figure(2);

worldmap([stn_latlon(1)-8 stn_latlon(1)+8],[stn_latlon(2)-8 stn_latlon(2)+8])
geoshow('landareas.shp','FaceColor','White')
title(sprintf('Earthquake Map | %s - %s, M > %0.1f, d < %.0f km',min(days_vec),max(days_vec),mag_min,dis_max));

for i=1:size(EQ_sel,1)
    EQ_lat=EQ_sel(i,8);
    EQ_lon=EQ_sel(i,9);
    EQ_radius=((EQ_sel(i,4)-4.5)/0.025);
    if EQ_sel(i,6)<=50
        colour='r';
    elseif (EQ_sel(i,6)>50) && (EQ_sel(i,6)<=200)
        colour='b';
    else EQ_sel(i,6)>200
         colour='g';
    end
    date_text=datetime(datevec(EQ_sel(i,1)),'Format','dd/MM/yyyy');
    circlem(EQ_lat,EQ_lon,EQ_radius,'edgecolor','none','facecolor',colour,'facealpha','0.2');
    if (EQ_sel(i,4)>=6.0 && EQ_sel(i,3)<=100) || EQ_sel(i,3)<=30 || EQ_sel(i,4)>=7.0
        textm(EQ_lat,EQ_lon,sprintf('%s(%0.1f)',date_text,EQ_sel(i,4)),'FontSize',8);
    end
end

plotm(stn_latlon(1),stn_latlon(2),'k^','MarkerFaceColor','k','MarkerSize',5);
textm(stn_latlon(1),stn_latlon(2)+0.2,stn);

for i=1:length(stn_MAGDAS)
    if stn_MAGDAS(i,2)>stn_latlon(1)-8 && stn_MAGDAS(i,2)<stn_latlon(1)+8 && stn_MAGDAS(i,3)>stn_latlon(2)-8 && stn_MAGDAS(i,3)<stn_latlon(2)+8
        stnnearby_latlon=stn_MAGDAS(i,2:3);
        stnnearby_name=stn_vec(i);
        plotm(stnnearby_latlon(1),stnnearby_latlon(2),'k^','MarkerSize',5);
        textm(stnnearby_latlon(1),stnnearby_latlon(2)+0.2,stnnearby_name);
    end
end

set(gcf, 'Position', get(0, 'Screensize'));

today=char(datetime('today','Format','dd-MM-yyyy'));
path1='D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\';
mkdir(fullfile(path1,today));

path2=strcat('D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\',today,'\');
figname=char(strcat(stn,year));
mapname=char(strcat(stn,year,' MAP'));
saveas(f1,fullfile(path2,figname),'tiff');
saveas(f2,fullfile(path2,mapname),'tiff');

