%% %Purpose
%To plot daily Z/H and geomagnetic indices temporal evolution
%% %Place and time inputs
stn='LWA';                         %Abbreviation of station name
date_start=[2015,09,18];           %Insert custom start and end dates
date_end=  [2015,09,23];           %Period spanning through 3 consecutive years is the maximum  
%% %Customization
mag_min=5.0;                     %Minimum magnitude of earthquakes to be considered
dis_max=500;                     %Maximum epicentral distance from the station

f_1=0.011;
f_2=0.013;

LT_start=22;
LT_end=02;
%% %Time period
datenum_start=datenum(date_start);
datenum_end=datenum(date_end);
if date_start(1)==date_end(1)
    year=string(date_start(1));
else
    year=sprintf('%s-%s',string(date_start(1)),string(date_end(1)));
end

load VARIABLES_WORLD

j=date_start(1);
for i=1:3
    year_vec(1,i)=j;
    if j==date_end(1)
        break
    end
    j=j+1;
end

days_num=datenum_end-datenum_start+1;
days_vec=datetime(datevec(datenum_start:1:datenum_end),'Format','dd/MM/yyyy');
%% %Loading files
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
%% %Station setting
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
    
stn_latlon=[stn_MAGDAS(stn_num,2:3)];
%% %Building earthquakes table
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
EQ_sel=NaN(1,7);
for j=1:size(EQ_table,1)
    if  EQ_table(j,10)>=mag_min
        EQ_latlon=[EQ_table(j,7:8)];
        EQ_dis=deg2km(distance('gc',stn_latlon,EQ_latlon));
        EQ_time=datenum(EQ_table(j,1:6));
        EQ_mag=EQ_table(j,10);
        EQ_Ks=(10^(0.75*EQ_mag))/(EQ_dis+100);
        EQ_depth=EQ_table(j,9);
       
        if EQ_dis<=dis_max
            EQ_sel(i,1)=EQ_time;
            EQ_sel(i,2)=EQ_dis;
            EQ_sel(i,3)=EQ_mag;
            EQ_sel(i,4)=EQ_Ks;
            EQ_sel(i,5)=EQ_depth;
            EQ_sel(i,6:7)=EQ_latlon;
            i=i+1;
        end
    end
end
%% %Geomagnetic data pre-processing
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

%Fill gaps
H_night=fillgaps(H_night,60);
D_night=fillgaps(D_night,60);
Z_night=fillgaps(Z_night,60);
%% %ap and Dst indices
hr_adjust=0;

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
%% %PSD calculation
n_win=3600;
n_ovrlp=0.5*n_win;
n_fft=n_win;
fs=1;
bin_1=round(((n_fft/2+1)/0.5)*f_1);
bin_2=round(((n_fft/2+1)/0.5)*f_2);

[H_wlch,f_wlch]=pwelch(H_night,hamming(n_win),n_ovrlp,n_fft,fs);
H_wlchf=H_wlch(bin_1:bin_2,:);
H_wlchdaymu=mean(H_wlchf,1,'omitnan');

[D_wlch,f_wlch]=pwelch(D_night,hamming(n_win),n_ovrlp,n_fft,fs);
D_wlchf=D_wlch(bin_1:bin_2,:);
D_wlchdaymu=mean(D_wlchf,1,'omitnan');

[Z_wlch,f_wlch]=pwelch(Z_night,hamming(n_win),n_ovrlp,n_fft,fs);
Z_wlchf=Z_wlch(bin_1:bin_2,:);
Z_wlchdaymu=mean(Z_wlchf,1,'omitnan');
%% %Normalisation
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
        H_monmu(1,k)=mean(H_wlchdaymu(i-j:i),'omitnan');
        D_monmu(1,k)=mean(D_wlchdaymu(i-j:i),'omitnan');
        Z_monmu(1,k)=mean(Z_wlchdaymu(i-j:i),'omitnan');
        
        H_monsig(1,k)=std(H_wlchdaymu(i-j:i),'omitnan');
        D_monsig(1,k)=std(D_wlchdaymu(i-j:i),'omitnan');
        Z_monsig(1,k)=std(Z_wlchdaymu(i-j:i),'omitnan');
        
        H_norm(1,i-j:i)=(H_wlchdaymu(1,i-j:i)-H_monmu(1,k))/H_monsig(1,k);
        D_norm(1,i-j:i)=(D_wlchdaymu(1,i-j:i)-D_monmu(1,k))/D_monsig(1,k);
        Z_norm(1,i-j:i)=(Z_wlchdaymu(1,i-j:i)-Z_monmu(1,k))/Z_monsig(1,k);
        
        k=k+1;
        j=0;
    end
    
end
%% %Power ratio
ZH=Z_norm./H_norm;

%mu+/-2sigma
ZH_mu=mean(ZH,'omitnan');
ZH_sig=std(ZH,'omitnan');
ZH_mup2sig=ones(1,days_num)*(ZH_mu+2*ZH_sig);
ZH_mum2sig=ones(1,days_num)*(ZH_mu-2*ZH_sig);
ZH_mm=movmean(ZH,5,'omitnan');
%% %Disturbed days and anomalies detection
dstb=NaN(1,days_num);
anom=NaN(1,days_num);
for i=1:days_num
    if ZH(1,i)>mean(ZH_mup2sig) || ZH(1,i)<mean(ZH_mum2sig)
        if ap(i,1)>50 || Dst(i,1)<-50
            dstb(1,i)=ZH(1,i);
            if i<days_num
                dstb(1,i+1)=ZH(1,i+1);
            end
            if i>1
                dstb(1,i-1)=ZH(1,i-1);
            end
        else
            anom(1,i)=ZH(1,i);
             if i<days_num
                anom(1,i+1)=ZH(1,i+1);
            end
            if i>1
                anom(1,i-1)=ZH(1,i-1);
            end
        end
    end
end
%% %Figure 1
title_sup=sprintf('%s station, %s - %s',station,min(days_vec),max(days_vec));

f1=figure(1);
stitle=suptitle(sprintf('%s',title_sup));
stitle_pos =get(stitle,'position');
stitle_pos=[stitle_pos(1) stitle_pos(2)+0.01 stitle_pos(3)];
set(stitle,'position',stitle_pos);

subplot(3,1,1)
plot(days_vec,ap,'k');
hold on
plot(days_vec,Dst,'r');
plot(days_vec,ones(length(days_vec),1)*(50),'k--');
plot(days_vec,ones(length(days_vec),1)*(-50),'r--');
hold off
ylabel('Dst               ap');
xlim([min(days_vec) max(days_vec)])
if max(ap)>50 || min(Dst)<-50
    ylim([min(Dst)-10 max(ap)+10])
else
    ylim([-60 60])
end
yyaxis right
set(gca, 'YColor', 'b')
for i=1:size(EQ_sel,1)
    j=datetime(datevec(EQ_sel(i,1)),'Format','dd/MM/yyyy HH:mm:ss');
    plot(j,EQ_sel(i,4),'bo','MarkerFaceColor','b','MarkerSize',5);
    if EQ_sel(i,4)==max(EQ_sel(:,4))
        label=sprintf('%s M%0.1f \t%.0fKM %.0fKM AWAY',string(datetime(j,'Format','dd/MM/yyyy')),EQ_sel(i,3),EQ_sel(i,5),EQ_sel(i,2));
        text(j,EQ_sel(i,4),label,'VerticalAlignment','top','HorizontalAlignment','left')
    end
    hold on
end
hold off
xlim([min(days_vec) max(days_vec)])
ylabel('K_{LS}');
title('Geomagnetic (left) & local seismicity (right) indices');

subplot(3,1,2);
plot(days_vec,ZH,'b',days_vec,ZH_mum2sig,'b--',days_vec,ZH_mup2sig,'b--');
hold on
plot(days_vec,anom,'r');
plot(days_vec,dstb,'g');
ylabel('Z/H');
legend(sprintf('\\Deltaf=%.3f-%.3f Hz | \\DeltaT=%.0f-%.0f LT',f_1,f_2,LT_start,LT_end),'location','northwest');
xlim([min(days_vec) max(days_vec)])
hold off
title('Z/H polarisation ratio');

subplot(3,1,3);
time_start=datenum(horzcat(date_start,[0,0,0]));
time_end=datenum(horzcat(date_end,[23,59,59]));
time_vec=datetime(datevec(time_start:1/86400:time_end),'Format','dd/MM/yyyy HH:mm:ss');
plot(time_vec,H_dn)
yyaxis right
plot(time_vec,Z_dn)
xlim([min(time_vec) max(time_vec)])
title('Raw H (left) & raw Z (right)');
%% %Figure 2
f2=figure(2);

worldmap([stn_latlon(1)-8 stn_latlon(1)+8],[stn_latlon(2)-8 stn_latlon(2)+8])
geoshow('landareas.shp','FaceColor','White')
title(sprintf('Earthquake Map | %s - %s, M > %0.1f, d < %.0f km',min(days_vec),max(days_vec),mag_min,dis_max));

for i=1:size(EQ_sel,1)
    EQ_lat=EQ_sel(i,6);
    EQ_lon=EQ_sel(i,7);
    EQ_radius=((EQ_sel(i,3)-4.5)/0.025);
    if EQ_sel(i,5)<=50
        colour='r';
    elseif (EQ_sel(i,5)>50) && (EQ_sel(i,5)<=200)
        colour='b';
    else EQ_sel(i,5)>200
         colour='g';
    end
    date_text=datetime(datevec(EQ_sel(i,1)),'Format','dd/MM/yyyy');
    circlem(EQ_lat,EQ_lon,EQ_radius,'edgecolor','none','facecolor',colour,'facealpha','0.2');
    if (EQ_sel(i,3)>=6.0 && EQ_sel(i,2)<=100) || EQ_sel(i,2)<=30 || EQ_sel(i,3)>=7.0
        textm(EQ_lat,EQ_lon,sprintf('%s(%0.1f)',date_text,EQ_sel(i,3)),'FontSize',8);
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
%% %Azimuthal angle of source
date_input=input('Insert date which contains anomalous signal in ddmmyyyy format or 0 to skip ')
if date_input==0
    sprintf('Analysis ends')
else
    date_in=datenum(date_input,'ddmmyyyy');
    date_idx=find(datenum_vec==date_in);
    
    mat_HD=[(H_night(:,date_idx))';(D_night(:,date_idx))'];
    mat_Z=(Z_night(:,date_idx));
    AB=(inv(mat_HD*mat_HD'))*(mat_HD*mat_Z);
    AB_amp=sqrt(AB(1)^2+AB(2)^2);
    AB_theta=atand(AB(2)/AB(1));
    
    x_az=10*AB_amp*cosd(AB_theta);
    y_az=10*AB_amp*sind(AB_theta);
    linem([stn_latlon(1),stn_latlon(1)+x_az],[stn_latlon(2),stn_latlon(2)+y_az],'LineWidth',3);
end
%% %Saving figures
MMM_yyyy=upper(string(datetime(mean(welchnum_vec),'Format','MMMyyyy')));

today=char(datetime('today','Format','dd-MM-yyyy'));
today=strcat(today,'\BMKG_6\');
path1='D:\OneDrive\Belajar\Dr. Adib Yusof\Pembelajaran\Jakarta Julai 2018\Figures\';
mkdir(fullfile(path1,today));

path2=strcat('D:\OneDrive\Belajar\Dr. Adib Yusof\Pembelajaran\Jakarta Julai 2018\Figures\',today);

fq_1=regexprep(string(f_1),'[.]','');
fq_2=regexprep(string(f_2),'[.]','');

figname1=char(strcat(stn,'-',MMM_yyyy,' (ZH)'));
saveas(f1,fullfile(path2,figname1),'tiff');
saveas(f1,fullfile(path2,figname1),'fig');
