%% %Purpose
%To plot raw, diff, spectrogram, FFT and bandpass filtered data for both Z
%and H components
%% %Place and time inputs
stn='LWA';                                      %Abbreviation of station name
date_start=[2015,09,18,    00,00,00];           %Insert custom start and end dates
date_end=  [2015,09,25,    23,59,59];           %Period spanning through 3 consecutive years is the maximum  
title2=('Precursor to earthquake');
ID=('');
%% %Customization
%Filter range for FFT
bpf1=0.023;
bpf2=0.025;

%Bandpass filter range for filtered data
bpf3=bpf1;
bpf4=bpf2;

mag_min=5.0;                     %Minimum magnitude of earthquakes to be considered
dis_max=300;                     %Maximum epicentral distance from the station
%% %Time vector creation
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

days_num=datenum_end-datenum_start+1;
time_vec=datetime(datevec(datenum_start:1/86400:datenum_end),'Format','dd/MM/yyyy HH:mm:ss');

%Load other variables/indices
load VARIABLES_WORLD
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
%% %Station and earthquake parameters
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

%Building selected earthquakes table

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

H_dn=fillgaps(H_dn,60);
D_dn=fillgaps(D_dn,60);
Z_dn=fillgaps(Z_dn,60);
%% %Acquiring ap and Dst indices 
j=1;
k=1;
for i=1:length(index_geomag)
    m=datenum(index_geomag(i,1),1,index_geomag(i,2));
    if m>=datenum_start && m<=datenum_end
       
        ap(k:k+3599,1)=index_geomag(i,6);
        Dst(k:k+3599,1)=index_geomag(i,5);
        k=k+3600;
    end
end
%% %Calcualtion of diff, FFT and bandpass filter
H_diff=diff(H_dn);
H_diff(length(time_vec))=NaN;
D_diff=diff(D_dn);
D_diff(length(time_vec))=NaN;
Z_diff=diff(Z_dn);
Z_diff(length(time_vec))=NaN;

L=length(H_dn);
n=2^nextpow2(L);
Fs=1;
ff = Fs*(0:(L/2))/L;

n_win=4096;
n_ovrlp=0.75*n_win;
n_fft=n_win;

Y=fft(filter1('bp',H_dn,'fc',[bpf1 bpf2]),n);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
H_fft=P1;

Y=fft(filter1('bp',Z_dn,'fc',[bpf1 bpf2]),n);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Z_fft=P1;

H_filt=filter1('bp',H_dn,'fc',[bpf3 bpf4]);
D_filt=filter1('bp',D_dn,'fc',[bpf3 bpf4]);
Z_filt=filter1('bp',Z_dn,'fc',[bpf3 bpf4]);
%% %Plotting figures

%Plotting seismic index chart
title1=sprintf('%s station, %s',station,datetime(mean(time_vec),'Format','MMM yyyy'));
title3=strcat(title1,' (',title2,')');

%Figure 1

f1=figure(1)
stitle=suptitle(sprintf('%s | H component',title3));
stitle_pos =get(stitle,'position');
stitle_pos=[stitle_pos(1) stitle_pos(2)+0.03 stitle_pos(3)];

subplot(3,2,1)
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

subplot(3,2,2)
plot(time_vec,H_dn);
ylabel('H_{raw}')
xlim([min(time_vec) max(time_vec)])
title('Raw H component');
% set(gca, 'XTickLabel',[]);

n_win=4096;
n_ovrlp=0.75*n_win;
n_fft=n_win;

subplot(3,2,3)
[sH,fH,tH]=spectrogram(Z_dn,'yaxis',hamming(n_win),n_ovrlp,n_fft,1);
spectrogram(H_dn,'yaxis',hamming(n_win),n_ovrlp,n_fft,1) 
caxis([0 30])
ylim([0 60])
colormap jet
colorbar('off')
colorbar('east')
title('Spectrogram of H');

subplot(3,2,4)
plot(time_vec,H_diff);
ylabel('H_{diff}')
ylim([-5 5])
xlim([min(time_vec) max(time_vec)]);
title('Diff(H)');


subplot(3,2,5)
plot(time_vec,H_filt,'r')
ylabel('H_{filtered}');
xlim([min(time_vec) max(time_vec)])
title(sprintf('Bandpass filtered of H (%0.3f - %.3f Hz)',bpf3,bpf4));

subplot(3,2,6)
plot(ff,H_fft)
xlabel('f (Hz)')
ylabel('H_{FFT}')
ylim([0 1.5*max(H_fft)])
xlim([0 0.1])
title(sprintf('FFT of H (filtered %0.3f - %.3f Hz)',bpf1,bpf2));
% title(sprintf('FFT of H (filtered >%0.3f Hz)',bpf1));

set(gcf, 'Position', get(0, 'Screensize'));

%Figure 2

f2=figure(2)
stitle=suptitle(sprintf('%s | Z component',title3));
stitle_pos =get(stitle,'position');
stitle_pos=[stitle_pos(1) stitle_pos(2)+0.03 stitle_pos(3)];

subplot(3,2,1)
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

subplot(3,2,2)
plot(time_vec,Z_dn);
ylabel('Z_{raw}')
xlim([min(time_vec) max(time_vec)])
title('Raw Z component');

subplot(3,2,3)
[sZ,fZ,tZ]=spectrogram(Z_dn,'yaxis',hamming(n_win),n_ovrlp,n_fft,1);
spectrogram(Z_dn,'yaxis',hamming(n_win),n_ovrlp,n_fft,1)
caxis([0 30])
ylim([0 60])
colormap jet
colorbar('off')
colorbar('east')
title('Spectrogram of Z');

subplot(3,2,4)
plot(time_vec,Z_diff);
ylabel('Z_{diff}')
ylim([-5 5])
xlim([min(time_vec) max(time_vec)])
title('Diff(Z)');

subplot(3,2,5)
plot(time_vec,Z_filt,'r')
ylabel('Z_{filtered}');
xlim([min(time_vec) max(time_vec)])
title(sprintf('Bandpass filtered of Z (%0.3f - %.3f Hz)',bpf3,bpf4));

subplot(3,2,6)
plot(ff,Z_fft)
xlabel('f (Hz)')
ylabel('Z_{FFT}')
xlim([0 0.1])
ylim([0 1.5*max(Z_fft)])
title(sprintf('FFT of Z (filtered %0.3f - %.3f Hz)',bpf1,bpf2));
% title(sprintf('FFT of Z (filtered >%0.3f Hz)',bpf1));

set(gcf, 'Position', get(0, 'Screensize'));
%--------------------------------------------------------------------------
%% %Saving files
MMM_yyyy=upper(string(datetime(mean(time_vec),'Format','MMMyyyy')));

today=char(datetime('today','Format','dd-MM-yyyy'));
today=strcat(today,'\BMKG_1\');
path1='D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\';
mkdir(fullfile(path1,today));

path2=strcat('D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\',today);
figname1=char(strcat(stn,'-',MMM_yyyy,'-',upper(title2),ID,'(H)'));
figname2=char(strcat(stn,'-',MMM_yyyy,'-',upper(title2),ID,'(Z)'));
saveas(f1,fullfile(path2,figname1),'tiff');
saveas(f2,fullfile(path2,figname2),'tiff');
saveas(f1,fullfile(path2,figname1),'fig');
saveas(f2,fullfile(path2,figname2),'fig');
%% %Azimuthal angle of source
mat_HD=[H_dn';D_dn'];
mat_Z=Z_dn;
AB=(inv(mat_HD*mat_HD'))*(mat_HD*mat_Z);

AB_amp=sqrt(AB(1)^2+AB(2)^2);
AB_theta=atand(AB(2)/AB(1));