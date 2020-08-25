%% %Purpose
%To determine a specific frequency that produces the most frequent maximum
%amplitude of PSD spectrum.
%% %Place and time inputs
stn='CEB';                                 %Abbreviation of station name
date_start=[2012,01,20,00,00,00];           %Insert custom start and end dates
date_end=  [2012,02,06,23,59,59];           %Period spanning through 3 consecutive years is the maximum  
%% %Customization
%Set magnitude % distance
mag_min=5.0;                     
dis_max=300;    

%Segmentation (in secs)
N_seg=3600;

%Range #1------------------------------------------------------------------
%Frequency range for PSD
f_1=0.007;
f_2=0.022;
%Bandpass filter range for filtered data
bpf1=f_1;
bpf2=f_2;

%Range #2------------------------------------------------------------------
%Frequency range for PSD
f_3=0.022;
f_4=0.100;
%Bandpass filter range for filtered data
bpf3=f_3;
bpf4=f_4;
%% %Load files and time calculation
%Loading files
load VARIABLES_WORLD

%Time period
datenum_start=datenum(date_start);
datenum_end=datenum(date_end);
if date_start(1)==date_end(1)
    year=string(date_start(1));
else
    year=sprintf('%s-%s',string(date_start(1)),string(date_end(1)));
end

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

%Creating time vector
time_vec=datetime(datevec(datenum_start:1/86400:datenum_end),'Format','dd/MM/yyyy HH:mm:ss');
timenum_vec=datenum_start:1/86400:datenum_end;
welchnum_vec=datenum_start:1/(86400/N_seg):datenum_end;
datenum_vec=datenum_start:1:datenum_end;
%% %Station setting 

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
%% %Building selected earthquakes table
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

j=1;
for i=1:length(H)
    if UT1m(1,i)==datenum_start 
        H_period(j:j+length(time_vec)-1,1)=H(i:i+length(time_vec)-1,1);
        D_period(j:j+length(time_vec)-1,1)=D(i:i+length(time_vec)-1,1);
        Z_period(j:j+length(time_vec)-1,1)=Z(i:i+length(time_vec)-1,1);     
        break
    end
end

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

%Fill gaps
H_dn=fillgaps(H_dn,60);
D_dn=fillgaps(D_dn,60);
Z_dn=fillgaps(Z_dn,60);
%% %Set nighttime data
% %Nighttime data in LT
% LT_length=LT_end-LT_start;
% if LT_length<0 
%     LT_length=LT_length+24;  
% end
% t_zone=timezone(stn_latlon(2),'degrees');
% t_start=LT_start+t_zone;      
% if t_start>24
%     t_start=t_start-24;
% end
% if t_start<0
%    t_start=t_start+24;
% end
% t_starts=t_start*3600;      
% t_int=LT_length*3600;
% 
% for i=1:days_num
%     
%     H_night(:,i)=H_dn(t_starts:t_starts+t_int-1);
%     D_night(:,i)=D_dn(t_starts:t_starts+t_int-1);
%     Z_night(:,i)=Z_dn(t_starts:t_starts+t_int-1);
%     t_starts=t_starts+86400;
% end
%% %PSD calculation

H_seg=reshape(H_dn,N_seg,[]);
D_seg=reshape(D_dn,N_seg,[]);
Z_seg=reshape(Z_dn,N_seg,[]);

n_win=1800;
n_ovrlp=0.5*n_win;
n_fft=n_win;
fs=1;
bin_1=round(((n_fft/2+1)/0.5)*f_1);
bin_2=round(((n_fft/2+1)/0.5)*f_2);
bin_3=round(((n_fft/2+1)/0.5)*f_3);
bin_4=round(((n_fft/2+1)/0.5)*f_4);

[H_wlch,f_wlch]=pwelch(H_seg,hamming(n_win),n_ovrlp,n_fft,fs);
H_wlch1=H_wlch(bin_1:bin_2,:);
H_wlch2=H_wlch(bin_3:bin_4,:);

[Z_wlch,f_wlch]=pwelch(Z_seg,hamming(n_win),n_ovrlp,n_fft,fs);
Z_wlch1=Z_wlch(bin_1:bin_2,:);
Z_wlch2=Z_wlch(bin_3:bin_4,:);
%% %Power ratio
ZH_wlch1=Z_wlch1./H_wlch1;
ZH_wlch2=Z_wlch2./H_wlch2;

frange1=(((f_2-f_1)/size(ZH_wlch1,1))*(1:size(ZH_wlch1,1)))+f_1;
frange2=(((f_4-f_3)/size(ZH_wlch2,1))*(1:size(ZH_wlch2,1)))+f_3;

[max_wlch1,max_i1]=max(ZH_wlch1);
[max_wlch2,max_i2]=max(ZH_wlch2);

max_f1=(((f_2-f_1)/size(ZH_wlch1,1))*max_i1)+f_1;
max_f2=(((f_4-f_3)/size(ZH_wlch2,1))*max_i2)+f_3;
%% %Number of occurence
j=1;
for i=1:length(frange1)
    max_occ1(1,i)=sum(max_f1(1,:)==frange1(1,i));
end

j=1;
for i=1:length(frange2)
    max_occ2(1,i)=sum(max_f2(1,:)==frange2(1,i));
end
%% %Acquiring ap and Dst indices
k=1;
for i=1:length(index_geomag)
    m=datenum(index_geomag(i,1),1,index_geomag(i,2));
    if m>=datenum_start && m<=datenum_end
        ap(k:k+3599,1)=index_geomag(i,6);
        Dst(k:k+3599,1)=index_geomag(i,5);
        k=k+3600;
    end
end
%% %Figure 1

title_sup=sprintf('%s station, %s',station,datetime(mean(time_vec),'Format','MMM yyyy'));

f1=figure(1);
stitle=suptitle(sprintf('%s',title_sup));
stitle_pos =get(stitle,'position');
stitle_pos=[stitle_pos(1) stitle_pos(2)+0.01 stitle_pos(3)];
set(stitle,'position',stitle_pos);

subplot(2,2,[1 2])
plot(time_vec,ap,'b');
hold on
plot(time_vec,Dst,'r');
plot(time_vec,ones(length(time_vec),1)*(50),'b--');
plot(time_vec,ones(length(time_vec),1)*(-50),'r--');
hold off
ylabel('Dst          ap');
xlim([min(time_vec) max(time_vec)])
if max(ap)>50 || min(Dst)<-50
    ylim([min(Dst)-10 max(ap)+10])
else
    ylim([-60 60])
end
yyaxis right
for i=1:size(EQ_sel,1)
    j=datetime(datevec(EQ_sel(i,1)),'Format','dd/MM/yyyy HH:mm:ss');
    plot(j,EQ_sel(i,4),'ko','MarkerFaceColor','k','MarkerSize',5);
    if EQ_sel(i,4)==max(EQ_sel(:,4))
        label=sprintf('%s M%0.1f %.0fKM',string(datetime(j,'Format','dd/MM/yyyy')),EQ_sel(i,3),EQ_sel(i,5));
        text(j,EQ_sel(i,4),label,'VerticalAlignment','top','HorizontalAlignment','left')
    end
    hold on
end
hold off
xlim([min(time_vec) max(time_vec)])
ylabel('K_{LS}');
title('Geomagnetic (left) & local seismicity (right) indices');

colpal=jet(size(max_wlch1,2));

subplot(2,2,3)
for i=1:size(max_wlch1,2)
    plot(max_f1(1,i),max_wlch1(1,i),'o','Color',colpal(i,:))
    hold on
end
hold off
ylabel('Z/H_{Max power}')
xlim([min(frange1) max(frange1)])
ylim([0 1.5*max(max_wlch1)])
xlabel('Frequency (Hz)')
title(sprintf('Z/H maximum power (%0.3f - %.3f Hz)',f_1,f_2));

subplot(2,2,4)
for i=1:size(max_wlch2,2)
    plot(max_f2(1,i),max_wlch2(1,i),'o','Color',colpal(i,:))
    hold on
end
hold off
ylabel('Z/H_{Max power}')
xlim([min(frange2) max(frange2)])
ylim([0 1.5*max(max_wlch2)])
xlabel('Frequency (Hz)')
title(sprintf('Z/H maximum power (%0.3f - %.3f Hz)',f_3,f_4));

set(gcf, 'Position', get(0, 'Screensize'));
%% %Figure 2
f2=figure(2);
stitle=suptitle(sprintf('%s',title_sup));
stitle_pos =get(stitle,'position');
stitle_pos=[stitle_pos(1) stitle_pos(2)+0.01 stitle_pos(3)];
set(stitle,'position',stitle_pos);

subplot(3,1,1)
plot(time_vec,ap,'b');
hold on
plot(time_vec,Dst,'r');
plot(time_vec,ones(length(time_vec),1)*(50),'b--');
plot(time_vec,ones(length(time_vec),1)*(-50),'r--');
hold off
ylabel('Dst          ap');
xlim([min(time_vec) max(time_vec)])
if max(ap)>50 || min(Dst)<-50
    ylim([min(Dst)-10 max(ap)+10])
else
    ylim([-60 60])
end
yyaxis right
for i=1:size(EQ_sel,1)
    j=datetime(datevec(EQ_sel(i,1)),'Format','dd/MM/yyyy HH:mm:ss');
    plot(j,EQ_sel(i,4),'ko','MarkerFaceColor','k','MarkerSize',5);
    if EQ_sel(i,4)==max(EQ_sel(:,4))
        label=sprintf('%s M%0.1f %.0fKM',string(datetime(j,'Format','dd/MM/yyyy')),EQ_sel(i,3),EQ_sel(i,5));
        text(j,EQ_sel(i,4),label,'VerticalAlignment','bottom','HorizontalAlignment','left')
    end
    hold on
end
hold off
xlim([min(time_vec) max(time_vec)])
ylabel('K_{LS}');
title('Geomagnetic (left) & local seismicity (right) indices');

subplot(3,1,2)
plot(time_vec,ZH_2);
ylabel('Power_{Z/H}')
xlim([min(time_vec) max(time_vec)])
ylim([min(ZH_2) 1.2*max(ZH_2)])
title(sprintf('Z/H power ratio (%0.3f - %.3f Hz)',f_3,f_4));
amp_max=max(ZH_2);
label=sprintf('%0.2f',amp_max);
for i=1:length(ZH_2)
    if ZH_2(i)==amp_max
        text(time_vec(i),amp_max,label,'VerticalAlignment','bottom','HorizontalAlignment','left')
        break
    end
end

subplot(3,1,3)
plot(time_vec,ZH2_filt,'r')
ylabel('Z/H_{filtered}');
xlim([min(time_vec) max(time_vec)])
title(sprintf('Bandpass filtered of Z/H (%0.3f - %.3f Hz)',bpf3,bpf4));

set(gcf, 'Position', get(0, 'Screensize'));
%--------------------------------------------------------------------------
%% %Saving figures
MMM_yyyy=upper(string(datetime(mean(time_vec),'Format','MMMyyyy')));

today=char(datetime('today','Format','dd-MM-yyyy'));
today=strcat(today,'\BMKG_3\');
path1='D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\';
mkdir(fullfile(path1,today));

path2=strcat('D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\',today);

fq_1=regexprep(string(f_1),'[.]','');
fq_2=regexprep(string(f_2),'[.]','');
fq_3=regexprep(string(f_3),'[.]','');
fq_4=regexprep(string(f_4),'[.]','');

figname1=char(strcat(stn,'-',MMM_yyyy,'-',fq_1,'-',fq_2,'(ZH)'));
saveas(f1,fullfile(path2,figname1),'tiff');
figname2=char(strcat(stn,'-',MMM_yyyy,'-',fq_3,'-',fq_4,'(ZH)'));
saveas(f2,fullfile(path2,figname2),'tiff');