%Place and time inputs-----------------------------------------------------

stn='ANC';                         %Abbreviation of station name
date_start=[2007,07,01];           %Insert custom start and end dates
date_end=  [2007,09,15];           %Period spanning through 3 consecutive years is the maximum

%Customization-------------------------------------------------------------

mag_min=5.0;                     %Minimum magnitude of earthquakes to be considered
dis_max=300;                     %Maximum epicentral distance from the station

f_A1=0.022;
f_A2=0.100;

LT_start=00;
LT_end=06;

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

%Nighttime data in LT
LT_length=LT_end-LT_start;
if LT_length<0 
    LT_length=LT_length+24;  
end
t_zone=timezone(stn_latlon(2),'degrees');
t_start=LT_start+t_zone;      
if t_start>24
    t_start=t_start-24;
end
if t_start<0
   t_start=t_start+24;
end
t_starts=t_start*3600;      
t_int=LT_length*3600;

for i=1:days_num
    
    H_night(:,i)=H_dn(t_starts:t_starts+t_int-1);
    D_night(:,i)=D_dn(t_starts:t_starts+t_int-1);
    Z_night(:,i)=Z_dn(t_starts:t_starts+t_int-1);
    t_starts=t_starts+86400;
    
end

%Getting ap and Dst indices
hr_adjust=1;

j=1;
k=1;
l=1;
for i=1:length(index_geomag)
    m=datenum(index_geomag(i,1),1,index_geomag(i,2));
    if m>=datenum_start && m<=datenum_end
        if index_geomag(i,3)==t_start-hr_adjust
            ap(k,1)=mean(index_geomag(i:i+LT_length+hr_adjust,6));
            Dst(k,1)=mean(index_geomag(i:i+LT_length+hr_adjust,5));
            k=k+1;
        end
        if index_geomag(i,3)==0
            apall(l,1)=mean(index_geomag(i:i+23,6));
            Dstall(l,1)=mean(index_geomag(i:i+23,5));
            l=l+1;
        end
        j=j+1;
        
        if j==days_num*24+1
            break;
        end
    end
end


%---------------------------------------------------------

Fs=1;
nfft=length(H_night);
bin_1=round(((nfft/2+1)/0.5)*f_A1);
bin_2=round(((nfft/2+1)/0.5)*f_A2);

H_dft=fft(H_night,nfft,1);
H_dft=H_dft(1:nfft/2+1,:);
H_psd=(1/(Fs*nfft))*abs(H_dft).^2;
H_psd(2:end-1)=2*H_psd(2:end-1);
H_psdf=H_psd(bin_1:bin_2,:);
H_psddaymu=mean(H_psdf,1,'omitnan');

D_dft=fft(D_night,nfft,1);
D_dft=D_dft(1:nfft/2+1,:);
D_psd=(1/(Fs*nfft))*abs(D_dft).^2;
D_psd(2:end-1)=2*D_psd(2:end-1);
D_psdf=D_psd(bin_1:bin_2,:);
D_psddaymu=mean(D_psdf,1,'omitnan');

Z_dft=fft(Z_night,nfft,1);
Z_dft=Z_dft(1:nfft/2+1,:);
Z_psd=(1/(Fs*nfft))*abs(Z_dft).^2;
Z_psd(2:end-1)=2*Z_psd(2:end-1);
Z_psdf=Z_psd(bin_1:bin_2,:);
Z_psddaymu=mean(Z_psdf,1,'omitnan');

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
        H_monmu(1,k)=mean(H_psddaymu(i-j:i),'omitnan');
        D_monmu(1,k)=mean(D_psddaymu(i-j:i),'omitnan');
        Z_monmu(1,k)=mean(Z_psddaymu(i-j:i),'omitnan');
        
        H_monsig(1,k)=std(H_psddaymu(i-j:i),'omitnan');
        D_monsig(1,k)=std(D_psddaymu(i-j:i),'omitnan');
        Z_monsig(1,k)=std(Z_psddaymu(i-j:i),'omitnan');
        
        H1_norm(1,i-j:i)=(H_psddaymu(1,i-j:i)-H_monmu(1,k))/H_monsig(1,k);
        D1_norm(1,i-j:i)=(D_psddaymu(1,i-j:i)-D_monmu(1,k))/D_monsig(1,k);
        Z1_norm(1,i-j:i)=(Z_psddaymu(1,i-j:i)-Z_monmu(1,k))/Z_monsig(1,k);
        
        k=k+1;
        j=0;
    end
    
end

%----------------------------------

ZG1=Z1_norm./H1_norm;
% ZG1=Z1_norm./D1_norm;

%mu+/-2sigma
ZG1_mu=mean(ZG1,'omitnan');
ZG1_sig=std(ZG1,'omitnan');
ZG1_mup2sig=ones(1,days_num)*(ZG1_mu+2*ZG1_sig);
ZG1_mum2sig=ones(1,days_num)*(ZG1_mu-2*ZG1_sig);
ZG1_mm=movmean(ZG1,5,'omitnan');

%------------------------------------
%Disturbed days and anomalies detection

dstb=NaN(1,days_num);
anom=NaN(1,days_num);
for i=1:days_num
    if ZG1(1,i)>mean(ZG1_mup2sig) || ZG1(1,i)<mean(ZG1_mum2sig)
        if ap(i,1)>50 || Dst(i,1)<-50
            dstb(1,i)=ZG1(1,i);
            if i<days_num
                dstb(1,i+1)=ZG1(1,i+1);
            end
            if i>1
                dstb(1,i-1)=ZG1(1,i-1);
            end
        else
            anom(1,i)=ZG1(1,i);
             if i<days_num
                anom(1,i+1)=ZG1(1,i+1);
            end
            if i>1
                anom(1,i-1)=ZG1(1,i-1);
            end
        end
    end
end

%----------------------------------------------

for i=1:days_num
    H_filt(:,i)=filter1('bp',H_night(:,i),'fc',[f_A1 f_A2],'fs',1);
    D_filt(:,i)=filter1('bp',D_night(:,i),'fc',[f_A1 f_A2],'fs',1);
    Z_filt(:,i)=filter1('bp',Z_night(:,i),'fc',[f_A1 f_A2],'fs',1);
end


H_fm=mean(H_filt,1);
D_fm=mean(D_filt,1);
Z_fm=mean(Z_filt,1);
ZG=Z_fm./(sqrt(H_fm.^2+D_fm.^2));

subplot(2,1,1)
plot(days_vec,ZG)
subplot(2,1,2)
plot(days_vec,ZG1)

for i=1:days_num
    H_filt(:,i)=filter1('bp',H_night(:,i),'fc',[f_A1 f_A2],'fs',1);
    D_filt(:,i)=filter1('bp',D_night(:,i),'fc',[f_A1 f_A2],'fs',1);
    H_hil(:,i)=hilbert(H_filt(:,i));
    D_hil(:,i)=hilbert(D_filt(:,i));
end

for i=1:size(H_hil,2)
    
    theta(:,i)=(atan(((2*abs(H_hil(:,i)).*abs(D_hil(:,i)))./((abs(D_hil(:,i)).^2)-(abs(H_hil(:,i)).^2))).*cos(angle(H_hil(:,i))-angle(D_hil(:,i)))))/2;
    
end

alpha=rad2deg(pi-theta);
alpha_mu=mean(alpha,1);

%---------------------------------------
%Plotting seismic index chart

f1=figure(1);

subplot(4,1,1)
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
ylabel('K_{LS}');
h=zeros(3,1);
h(1) = plot(days_vec,NaN(1,days_num),'ro');
h(2) = plot(days_vec,NaN(1,days_num),'b^');
h(3) = plot(days_vec,NaN(1,days_num),'gx');
lgd=legend(h,'<50km','50-200km','>200km','location','northwest');
title(lgd,'Depth');
title(sprintf('Z/G polarisation ratio | %s station, %s - %s',station,min(days_vec),max(days_vec)));
hold off

%Plotting polarisation ratio

subplot(4,1,2);
plot(days_vec,ZG1,'b',days_vec,ZG1_mum2sig,'b--',days_vec,ZG1_mup2sig,'b--');
hold on
plot(days_vec,anom,'r');
plot(days_vec,dstb,'g');
ylabel('Z/G');
legend(sprintf('\\Deltaf=%.3f-%.3f Hz | \\DeltaT=%.0f-%.0f LT',f_A1,f_A2,LT_start,LT_end),'location','northwest');
xlim([min(days_vec) max(days_vec)])
set(gca, 'XTickLabel',[]);
hold off

%Plotting ap and Dst indices
subplot(4,1,3)
plot(days_vec,ap,'b');
hold on
plot(days_vec,Dst,'r');
plot(days_vec,ones(days_num,1)*(50),'b--');
plot(days_vec,ones(days_num,1)*(-50),'r--');
hold off
set(gca, 'XTickLabel',[]);
ylabel('Dst     ap');
legend('ap','Dst','location','northwest');
xlim([min(days_vec) max(days_vec)])
ylim([min(Dst)-10 max(ap)+10])

%Plotting phase angles
subplot(4,1,4)
plot(days_vec,alpha_mu)
xlabel('Date')
ylabel('Azimuthal angle (\theta)')

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
% 
% today=char(datetime('today','Format','dd-MM-yyyy'));
% today=strcat(today,'\Adib_8\');
% path1='D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\';
% mkdir(fullfile(path1,today));
% 
% path2=strcat('D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\',today);
% figname=char(strcat(stn,year));
% mapname=char(strcat(stn,year,' MAP'));
% saveas(f1,fullfile(path2,figname),'tiff');
% saveas(f2,fullfile(path2,mapname),'tiff');
% 
