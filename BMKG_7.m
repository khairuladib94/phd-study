%% %Purpose
%To plot daily Z/H and geomagnetic indices temporal evolution with
%customizable segmentation period (hourly as default) and nighttime period labels
%OPTIONAL:
%-Calculate azimuthal angle and amplitude of anomalies
%-Inspect raw data
%% %Place and time inputs
stn='PTK';                         %Abbreviation of station name
date_start=[2013,05,01];           %Insert custom start and end dates
date_end=  [2013,05,31];           %Period spanning through 3 consecutive years is the maximum  
%% %Customization
mag_min=5.0;                     %Minimum magnitude of earthquakes to be considered
dis_max=500;                     %Maximum epicentral distance from the station

f_1=0.010;
f_2=0.012;

N_seg=3600;

LT_start=0;
LT_end=6;
%% %Time period

date_start1=datenum(horzcat(date_start,[0,0,0]));
date_end1=datenum(horzcat(date_end,[23,59,59]));
days_vec=datetime(datevec(date_start1:1/(86400/N_seg):date_end1),'Format','dd/MM/yyyy HH');
datenum_start=datenum(date_start1);
datenum_end=datenum(date_end1);
days_num=length(days_vec);

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
EQ_table=NaN(1,11);
j=1;
for i=1:length(EQ_WORLD)
    if datenum(EQ_WORLD(i,1:6))>=datenum_start && datenum(EQ_WORLD(i,1:6))<=datenum_end
        EQ_table(j,:)=EQ_WORLD(i,:);
        j=j+1;
    end
    if datenum(EQ_WORLD(i,1:6))>datenum_end
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
day_end=floor(datenum_end)-data_start+1;
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

% Fill gaps
H_dn=fillgaps(H_dn,3600,70);
D_dn=fillgaps(D_dn,3600,70);
Z_dn=fillgaps(Z_dn,3600,70);

H_seg=reshape(H_dn,N_seg,[]);
D_seg=reshape(D_dn,N_seg,[]);
Z_seg=reshape(Z_dn,N_seg,[]);
%% %Nighttime data in LT

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
t_end=t_start+LT_length;
%% %ap and Dst indices
k=1;
for i=1:length(index_geomag)
    m=datenum(index_geomag(i,1),1,index_geomag(i,2));
    if m>=datenum_start && m<=datenum_end
        ap(k,1)=index_geomag(i,6);
        Dst(k,1)=index_geomag(i,5);
        k=k+1;
    end
end
%% %PSD calculation
n_win=1800;
n_ovrlp=0.5*n_win;
n_fft=n_win;
fs=1;
bin_1=round(((n_fft/2+1)/0.5)*f_1);
bin_2=round(((n_fft/2+1)/0.5)*f_2);

[H_wlch,f_wlch]=pwelch(H_seg,hamming(n_win),n_ovrlp,n_fft,fs);
% H_wlch=10*log10(H_wlch);
H_wlchf=H_wlch(bin_1:bin_2,:);
H_wlchmu=mean(H_wlchf,1,'omitnan');

[D_wlch,f_wlch]=pwelch(D_seg,hamming(n_win),n_ovrlp,n_fft,fs);
% D_wlch=10*log10(D_wlch);
D_wlchf=D_wlch(bin_1:bin_2,:);
D_wlchmu=mean(D_wlchf,1,'omitnan');

[Z_wlch,f_wlch]=pwelch(Z_seg,hamming(n_win),n_ovrlp,n_fft,fs);
% Z_wlch=10*log10(Z_wlch);
Z_wlchf=Z_wlch(bin_1:bin_2,:);
Z_wlchmu=mean(Z_wlchf,1,'omitnan');

H_wlchmu(H_wlchmu==0)=NaN;
D_wlchmu(D_wlchmu==0)=NaN;
Z_wlchmu(Z_wlchmu==0)=NaN; 
%% %Power ratio
ZH=detrend(Z_wlchmu./H_wlchmu);

%mu+/-2sigma
ZH_mu=mean(ZH,'omitnan');
ZH_sig=std(ZH,'omitnan');
ZH_mup2sig=ones(1,days_num)*(ZH_mu+2*ZH_sig);
ZH_mum2sig=ones(1,days_num)*(ZH_mu-2*ZH_sig);
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
%% %Azimuthal angle calculation

for i=1:size(H_seg,2)
    mat_HD=[H_seg(:,i)';D_seg(:,i)'];
    mat_Z=(Z_seg(:,i));
    AB=(inv(mat_HD*mat_HD'))*(mat_HD*mat_Z);
    AB_amp(1,i)=sqrt(AB(1)^2+AB(2)^2);
    AB_theta(1,i)=atand(AB(2)/AB(1));
end
%% %Figure 1
title_sup=sprintf('%s station, %s - %s',station,datetime(datevec(datenum(date_start)),'Format','dd/MM/yyyy'),datetime(datevec(datenum(date_end)),'Format','dd/MM/yyyy'));

% f1=figure(1);
figure
stitle=suptitle(sprintf('%s',title_sup));
stitle_pos =get(stitle,'position');
stitle_pos=[stitle_pos(1) stitle_pos(2)+0.01 stitle_pos(3)];
set(stitle,'position',stitle_pos);

subplot(5,1,1)
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

subplot(5,1,2)
spectrogram(H_dn,'yaxis',hamming(n_win),n_ovrlp,n_fft,1)
caxis([0 10])
ylim([15 60])
colormap jet
colorbar('off')
colorbar('east')
title('Spectrogram of H');

subplot(5,1,3)
spectrogram(Z_dn,'yaxis',hamming(n_win),n_ovrlp,n_fft,1)
caxis([0 10])
ylim([15 60])
colormap jet
colorbar('off')
colorbar('east')
title('Spectrogram of Z');

subplot(5,1,4);
plot(days_vec,ZH,'b',days_vec,ZH_mum2sig,'k--',days_vec,ZH_mup2sig,'k--');
hold on
for i=0:size(days_vec,1)/(86400/N_seg)
    xpatch=[t_start/24+i t_end/24+i t_end/24+i t_start/24+i];
    ypatch=[-1e+5 -1e+5 1e+5 1e+5];
    patch(xpatch,ypatch,'k','FaceAlpha',0.05,'EdgeAlpha',0);
end
plot(days_vec,ZH,'b',days_vec,ZH_mum2sig,'k--',days_vec,ZH_mup2sig,'k--');
plot(days_vec,anom,'r');
plot(days_vec,dstb,'g');
plot(days_vec,movmean(ZH,10),'m','LineWidth',2)
hold off
ylabel('Z/H power ratio');
legend(sprintf('\\Deltaf=%.3f-%.3f Hz',f_1,f_2),'location','northwest');
xlim([min(days_vec) max(days_vec)])
ylim([min(ZH)+mean(ZH_mum2sig) max(ZH)+mean(ZH_mup2sig)])
title('Z/H polarisation ratio');

subplot(5,1,5)
plot(days_vec,movmean(AB_theta,12))
ylabel('Azimuthal angle (\theta)');
yyaxis right
plot(days_vec,AB_amp)
ylabel('Amplitude');
xlim([min(days_vec) max(days_vec)])

set(gcf, 'Position', get(0, 'Screensize'));
%% %Figure 2
f2=figure(2);

worldmap([stn_latlon(1)-8 stn_latlon(1)+8],[stn_latlon(2)-8 stn_latlon(2)+8])
geoshow('landareas.shp','FaceColor','White')
title(sprintf('Earthquake Map | %s - %s, M > %0.1f, d < %.0f km',datetime(datevec(datenum(date_start)),'Format','dd/MM/yyyy'),datetime(datevec(datenum(date_end)),'Format','dd/MM/yyyy'),mag_min,dis_max));

if ~isnan(mean(EQ_sel))
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

% %% %Azimuthal angle of source
% date_input=input('Insert date which contains anomalous signal in ddmmyyyy format or 0 to skip ')
% if date_input==0
%     sprintf('Analysis ends')
% else
%     date_in=datenum(date_input,'ddmmyyyy');
%     date_idx=find(datenum_vec==date_in);
%     
%     mat_HD=[(H_night(:,date_idx))';(D_night(:,date_idx))'];
%     mat_Z=(Z_night(:,date_idx));
%     AB=(inv(mat_HD*mat_HD'))*(mat_HD*mat_Z);
%     AB_amp=sqrt(AB(1)^2+AB(2)^2);
%     AB_theta=atand(AB(2)/AB(1));
%     
%     x_az=10*AB_amp*cosd(AB_theta);
%     y_az=10*AB_amp*sind(AB_theta);
%     linem([stn_latlon(1),stn_latlon(1)+x_az],[stn_latlon(2),stn_latlon(2)+y_az],'LineWidth',3);
% end
%% %Saving figures
MMM_yyyy=upper(string(datetime(mean(days_vec),'Format','MMMyyyy')));

today=char(datetime('today','Format','dd-MM-yyyy'));
today=strcat(today,'\BMKG_7\');
path1='D:\OneDrive\Belajar\Dr. Adib Yusof\Pembelajaran\Jakarta Julai 2018\Figures\';
mkdir(fullfile(path1,today));

path2=strcat('D:\OneDrive\Belajar\Dr. Adib Yusof\Pembelajaran\Jakarta Julai 2018\Figures\',today);

fq_1=regexprep(string(f_1),'[.]','');
fq_2=regexprep(string(f_2),'[.]','');
figname1=char(strcat(stn,MMM_yyyy,fq_1,'-',fq_2,'(ZH)'));
figname2=char(strcat(stn,MMM_yyyy,'(MAP)'));

saveas(f1,fullfile(path2,figname1),'tiff');
saveas(f1,fullfile(path2,figname1),'fig');

saveas(f2,fullfile(path2,figname2),'tiff');
saveas(f2,fullfile(path2,figname2),'fig');
%% %Inspect data
inspect=input('Inspect data? ')
if inspect==0
    sprintf('Analysis ends')
elseif inspect==1
    datenum_instart=date_start1;
    datenum_insend=date_end1;
else
    date_instart=inspect;
    inspect2=input('End date? ')
    if inspect2==1
        date_insend=date_instart;
        datenum_instart=datenum(date_instart,'ddmmyyyy');
        datenum_insend=datenum(horzcat(date_instart,'235959'),'ddmmyyyyhhMMss');
    else
        date_insend=inspect2;
        datenum_instart=datenum(date_instart,'ddmmyyyy');
        datenum_insend=datenum(horzcat(date_insend,'235959'),'ddmmyyyyhhMMss');
    end

end

if inspect~=0
    
instart=datenum_instart-data_start+1;
insend=floor(datenum_insend)-data_start+1;
insdatenum_vec=datenum_instart:datenum_insend;
insdays_vec=datetime(datevec(datenum_instart:1/86400:datenum_insend),'Format','dd/MM/yyyy HH:mm:ss');

H_ins=H(86400*instart-86399:86400*insend);
Z_ins=Z(86400*instart-86399:86400*insend);

%Removing noise/outlier
H_ins=medfilt1(H_ins,'omitnan');
Z_ins=medfilt1(Z_ins,'omitnan');

for i=1:5
    
    H_sig=std(H_ins,'omitnan');
    H_mu=mean(H_ins,'omitnan');
    Z_sig=std(Z_ins,'omitnan');
    Z_mu=mean(Z_ins,'omitnan');
    
    for j=1:length(H_ins)
 
        if H_ins(j)>H_mu+5*H_sig||H_ins(j)<H_mu-5*H_sig
            H_ins(j)=NaN;
        end
        if Z_ins(j)>Z_mu+5*Z_sig||Z_ins(j)<Z_mu-5*Z_sig
            Z_ins(j)=NaN;
        end
        
    end
    
end

H_diff=diff(H_ins);
H_diff(length(insdays_vec))=NaN;
Z_diff=diff(Z_ins);
Z_diff(length(insdays_vec))=NaN;

n_win=1800;
n_ovrlp=0.5*n_win;
n_fft=n_win;
fs=1;

title_sup=sprintf('Data inspection at %s station, %s - %s',station,datetime(datevec(datenum_instart),'Format','dd/MM/yyyy'),datetime(datevec(floor(datenum_insend)),'Format','dd/MM/yyyy'));

f3=figure(3);
stitle=suptitle(sprintf('%s',title_sup));
stitle_pos =get(stitle,'position');
stitle_pos=[stitle_pos(1) stitle_pos(2)+0.01 stitle_pos(3)];
set(stitle,'position',stitle_pos);

subplot(3,2,1)
plot(insdays_vec,H_ins);
ylabel('H_{raw}')
xlim([min(insdays_vec) max(insdays_vec)])
title('Raw H component');

subplot(3,2,2)
plot(insdays_vec,Z_ins);
ylabel('Z_{raw}')
xlim([min(insdays_vec) max(insdays_vec)])
title('Raw Z component');

subplot(3,2,3)
plot(insdays_vec,H_diff);
ylabel('H_{diff}')
ylim([-3 3])
xlim([min(insdays_vec) max(insdays_vec)])
title('Diff(H)');

subplot(3,2,4)
plot(insdays_vec,Z_diff);
ylabel('Z_{diff}')
ylim([-3 3])
xlim([min(insdays_vec) max(insdays_vec)])
title('Diff(Z)');

subplot(3,2,5)
spectrogram(H_ins,'yaxis',hamming(n_win),n_ovrlp,n_fft,1)
caxis([0 30])
ylim([0 60])
colormap jet
colorbar('off')
colorbar('east')
title('Spectrogram of H');

subplot(3,2,6)
spectrogram(Z_ins,'yaxis',hamming(n_win),n_ovrlp,n_fft,1)
caxis([0 30])
ylim([0 60])
colormap jet
colorbar('off')
colorbar('east')
title('Spectrogram of Z');
    
set(gcf, 'Position', get(0, 'Screensize'));

end





 try3=8.476*exp(-111.1*try2)
 subplot(2,1,1)
 plot(try2,try1)
 hold on
 plot(try2,try3,'r')
 hold off
 
 subplot(2,1,2)
 plot(try2,try4)
 

 
