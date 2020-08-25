%% Purpose
% Using polarization ellipse methodas in Ohta et. al (2013)
% Search for ULF precursor in frequency bands. Frequency might have been
% changed to higher ones (>0.1 Hz); check f_rangeint.
% Uses only ONE local nighttime period
% Include phase shift for determinancy (deltaS) calculation
% Altered from Adib_105, where the segmentation of PSD calculation are
% strictly based on Schekotov 2008
%% Load initial variable
load VARIABLES_WORLD
open EQ_intrst
%% Place and date
EQ_num=2;          %Select earthquake of interest as listed in EQ_intrst #

stn0=EQ_intrst{EQ_num,1}; stn=stn0(1:3);
EQ_date0=EQ_intrst{EQ_num,3}; EQ_date=datenum(str2num(strcat('20',EQ_date0(7:8))),str2num(EQ_date0(4:5)),str2num(EQ_date0(1:2)));

date_start0=EQ_date-45;
date_end0=EQ_date+15;

date_start=datevec(date_start0);
date_end=datevec(date_end0);
%% Customization
mag_min=5.0;                     %Minimum magnitude of earthquakes to be considered
dis_max=300;                     %Maximum epicentral distance from the station
depth_max=200;

LT_start=01;
LT_end=05;

f_rangeint=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10];

k_Thres=5; N_medfil=1;  k_repeat=5;     %Noise removal
k_alpha=4;                              %Minimal SNR coefficient
%% Time period
date_start1=datenum(horzcat(date_start(1:3),[0,0,0]));
date_end1=datenum(horzcat(date_end(1:3),[23,59,59]));
datenum_start=datenum(date_start);
datenum_end=datenum(date_end);
days_num=datenum_end-datenum_start+1;
days_vec=datetime(datevec(datenum_start:1:datenum_end),'Format','dd/MM/yyyy');

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

%% Loading files
filepathname=strcat('E:\Study\MAGDAS DATA\',stn,'\',stn,string(year_vec(:)),'S.mat');

if numel(year_vec)==1
    matname=strcat(stn,string(year_vec(1)),'S');
    load(matname);
end

if numel(year_vec)==2 
    matname(1)=strcat(stn,string(year_vec(1)),'S');
    matname(2)=strcat(stn,string(year_vec(2)),'S');
    if exist(filepathname(1),'file')==2 && exist(filepathname(2),'file')==2
        A=load(matname(1));
        B=load(matname(2));
    elseif exist(filepathname(1),'file')==2 && exist(filepathname(2),'file')==0
        A=load(matname(1));
        B.UT1m=datenum(year_vec(2),01,01,00,00,00):1/86400:datenum(year_vec(2),12,31,23,59,59);
        B.H=NaN(numel(B.UT1m),1); B.D=NaN(numel(B.UT1m),1); B.Z=NaN(numel(B.UT1m),1); B.F=NaN(numel(B.UT1m),1);
    elseif exist(filepathname(1),'file')==0 && exist(filepathname(2),'file')==2
        B=load(matname(2));
        A.UT1m=datenum(year_vec(1),01,01,00,00,00):1/86400:datenum(year_vec(1),12,31,23,59,59);
        A.H=NaN(numel(A.UT1m),1); A.D=NaN(numel(A.UT1m),1); A.Z=NaN(numel(A.UT1m),1); A.F=NaN(numel(A.UT1m),1);
    end
    H=vertcat(A.H,B.H);
    D=vertcat(A.D,B.D);
    Z=vertcat(A.Z,B.Z);
    UT1m=horzcat(A.UT1m,B.UT1m);
end
%% Station setting
for i=1:length(stn_MAGDAS)
    if strcmp(stn,stn_vec(i))
        stn_num=i;
        station=string(station_vec(i,1));
        region=string(station_vec(i,2));
        break;
    end
end

stn_latlon=[stn_MAGDAS(stn_num,2:3)];
%% Building earthquakes table
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
EQ_sel=NaN(1,11);
for j=1:size(EQ_table,1)
    if  EQ_table(j,10)>=mag_min
        EQ_latlon=[EQ_table(j,7:8)];
        EQ_dis=deg2km(distance('gc',stn_latlon,EQ_latlon));
        EQ_time=datenum(EQ_table(j,1:6));
        EQ_mag=EQ_table(j,10);
        EQ_Ks=(10^(0.75*EQ_mag))/(EQ_dis+100);
        EQ_depth=EQ_table(j,9);
        EQ_angle=azimuth('gc',stn_latlon(1:2),EQ_latlon(1:2),'degree');
        
        if EQ_dis<=dis_max && EQ_depth<=depth_max
            EQ_sel(i,1)=EQ_time;
            EQ_sel(i,2)=EQ_dis;
            EQ_sel(i,3)=EQ_mag;
            EQ_sel(i,4)=EQ_Ks;
            EQ_sel(i,5)=EQ_depth;
            EQ_sel(i,6:7)=EQ_latlon;
            EQ_sel(i,8)=EQ_angle;
            i=i+1;
        end
    end
end
colpal=jet(size(EQ_sel,1));
EQ_sel(:,9:11)=colpal;
%% Geomagnetic data pre-processing
data_start=round(min(UT1m));
data_end=floor(max(UT1m));

day_start=datenum_start-data_start+1;
day_end=floor(datenum_end)-data_start+1;
datenum_vec=datenum_start:datenum_end;
UT1m_sel=UT1m(86400*day_start-86399:86400*day_end);

H_dn=H(86400*day_start-86399:86400*day_end);
D_dn=D(86400*day_start-86399:86400*day_end);
Z_dn=Z(86400*day_start-86399:86400*day_end);
G_dn=sqrt(H_dn.^2+D_dn.^2);

H_ori=H(86400*day_start-86399:86400*day_end);
D_ori=D(86400*day_start-86399:86400*day_end);
Z_ori=Z(86400*day_start-86399:86400*day_end);
G_ori=sqrt(H_ori.^2+D_ori.^2);

%Removing noise/outlier

H_dn=medfilt1(H_dn,N_medfil,'omitnan','truncate');
D_dn=medfilt1(D_dn,N_medfil,'omitnan','truncate');
Z_dn=medfilt1(Z_dn,N_medfil,'omitnan','truncate');
G_dn=medfilt1(G_dn,N_medfil,'omitnan','truncate');

for i=1:k_repeat
    
    H_sig=std(H_dn,'omitnan');
    H_mu=mean(H_dn,'omitnan');
    D_sig=std(D_dn,'omitnan');
    D_mu=mean(D_dn,'omitnan');
    Z_sig=std(Z_dn,'omitnan');
    Z_mu=mean(Z_dn,'omitnan');
    G_sig=std(G_dn,'omitnan');
    G_mu=mean(G_dn,'omitnan');
    
    for j=1:length(H_dn)
        
        if H_dn(j)>H_mu+k_Thres*H_sig||H_dn(j)<H_mu-k_Thres*H_sig
            H_dn(j)=NaN;
        end
        if D_dn(j)>D_mu+k_Thres*D_sig||D_dn(j)<D_mu-k_Thres*D_sig
            D_dn(j)=NaN;
        end
        if Z_dn(j)>Z_mu+k_Thres*Z_sig||Z_dn(j)<Z_mu-k_Thres*Z_sig
            Z_dn(j)=NaN;
        end
        if G_dn(j)>G_mu+k_Thres*G_sig||G_dn(j)<G_mu-k_Thres*G_sig
            G_dn(j)=NaN;
        end
    end
    
end
%% Nighttime data in LT and gap-filling

%Nighttime segmentation

LT_length=LT_end-LT_start;
if LT_length<0 LT_length=LT_length+24; end
t_zone=timezone(stn_latlon(2),'degrees');
t_start=LT_start+t_zone;
if t_start>24 t_start=t_start-24; end
if t_start<0 t_start=t_start+24;end
t_starts=t_start*3600;
t_int=LT_length*3600;

for i=1:days_num
    H_night(:,i)=H_dn(t_starts+1:t_starts+t_int);
    D_night(:,i)=D_dn(t_starts+1:t_starts+t_int);
    Z_night(:,i)=Z_dn(t_starts+1:t_starts+t_int);
    G_night(:,i)=G_dn(t_starts+1:t_starts+t_int);
    UT1m_night(:,i)=UT1m_sel(t_starts+1:t_starts+t_int);
    t_starts=t_starts+86400;
end

%Gap-filling
for i=1:days_num
    H_measure=regionprops(logical(isnan(H_night(:,i))), 'Area');
    H_length=[H_measure.Area]; if isempty(H_length) H_length=1; end
    if max(H_length)<0.01*size(H_night,1) && any(isnan(H_night(:,i)))
        H_night(:,i)=fillgaps(H_night(:,i));
    end
    
    D_measure=regionprops(logical(isnan(D_night(:,i))), 'Area');
    D_length=[D_measure.Area]; if isempty(D_length) D_length=1; end
    if max(D_length)<0.01*size(D_night,1) && any(isnan(D_night(:,i)))
        D_night(:,i)=fillgaps(D_night(:,i));
    end
    
    Z_measure=regionprops(logical(isnan(Z_night(:,i))), 'Area');
    Z_length=[Z_measure.Area]; if isempty(Z_length) Z_length=1; end
    if max(Z_length)<0.01*size(Z_night,1) && any(isnan(Z_night(:,i)))
        Z_night(:,i)=fillgaps(Z_night(:,i));
    end
    
    G_measure=regionprops(logical(isnan(G_night(:,i))), 'Area');
    G_length=[G_measure.Area]; if isempty(G_length) G_length=1; end
    if max(G_length)<0.01*size(G_night,1) && any(isnan(G_night(:,i)))
        G_night(:,i)=fillgaps(G_night(:,i));
    end
 
end

%% ap and Dst indices
geomag_start0=find(datenum(index_geomag(:,1),1,index_geomag(:,2))==datenum_start);
geomag_start=geomag_start0(1);
geomag_end0=find(datenum(index_geomag(:,1),1,index_geomag(:,2))==datenum_end);
geomag_end=geomag_end0(end);
ap=index_geomag(geomag_start:geomag_end,6);
Dst=index_geomag(geomag_start:geomag_end,5);
hourly_vec=datetime(datevec(datenum_start:1/24:datenum_end+1),'Format','dd/MM/yyyy HH');
hourly_vec(end)=[];
%% Precursor and direction

alpha1=NaN(size(H_night,1),days_num);
alpha2=NaN(size(H_night,1),days_num);

for i=1:days_num
   
    %% DeltaSF (rotated coordinate)
    n_length=size(H_night,1);
    n_win=1024;
    n_ovrlp=0;
    n_fft=n_win;
    fs=1;
    
    % Filtering
    H_pass0=filter1('bp',H_night(:,i),'fc',[f_rangeint(1) f_rangeint(end)],'fs',fs,'order',4);
    D_pass0=filter1('bp',D_night(:,i),'fc',[f_rangeint(1) f_rangeint(end)],'fs',fs,'order',4);
    
    % Exclude excessive points
    excl_pts=rem(n_length,n_win);
    H_pass=H_pass0(1:end-excl_pts);  D_pass=D_pass0(1:end-excl_pts);

    % Power spectral density
    [Phh,fff]=pwelch(H_pass,hamming(n_win),n_ovrlp,n_fft,fs);
    [Pdd,fff]=pwelch(D_pass,hamming(n_win),n_ovrlp,n_fft,fs);
    
    [Phd,fff]=cpsd(H_pass,D_pass,hamming(n_win),n_ovrlp,n_fft,fs);
    [Pdh,fff]=cpsd(D_pass,H_pass,hamming(n_win),n_ovrlp,n_fft,fs);
     
    % ULF extraction  
    bin_floor=round(length(fff)/0.5*f_rangeint(1)); 
    bin_ceil=round(length(fff)/0.5*f_rangeint(end)); 
    
    Phh1=Phh(bin_floor:bin_ceil,:);
    Pdd1=Pdd(bin_floor:bin_ceil,:);
    
    Phd1=Phd(bin_floor:bin_ceil,:);
    Pdh1=Pdh(bin_floor:bin_ceil,:);
    
    fff1=linspace(f_rangeint(1),f_rangeint(end),length(Phh1));
    
    % ULF averaging
    for b=1:numel(f_rangeint) f_rangeidx(b)=find(round(fff1,4)==f_rangeint(b)); end
    
    for a=1:numel(f_rangeint)-1
        Phh_band(:,a)=Phh1(f_rangeidx(a):f_rangeidx(a+1)-1);
        Pdd_band(:,a)=Pdd1(f_rangeidx(a):f_rangeidx(a+1)-1);
        
        Phd_band(:,a)=Phd1(f_rangeidx(a):f_rangeidx(a+1)-1);
        Pdh_band(:,a)=Pdh1(f_rangeidx(a):f_rangeidx(a+1)-1);
        
        fff_band(:,a)=fff1(f_rangeidx(a):f_rangeidx(a+1)-1);
    end
    
    Phh_mean=mean(Phh_band,1);
    Pdd_mean=mean(Pdd_band,1);
    
    Phd_mean=mean(Phd_band,1);
    Pdh_mean=mean(Pdh_band,1);
    
    % DeltaSF calculations
    beta=0.5*asin(imag(Pdh_mean-Phd_mean) ./ ((Phh_mean-Pdd_mean).^2 + 4*Phh_mean.*Pdd_mean ).^0.5);
    deltaS(:,i)=real(((Phh_mean ./ Pdd_mean) - 1) ./ rms(tan(beta)));
    
    del_psi=10; psi=deg2rad(0:del_psi:180-del_psi);
    for f=1:numel(Phh_mean)
        Ptt=Phh_mean(f).*(cos(psi)).^2 + Pdd_mean(f).*(sin(psi)).^2 + real(Phd_mean(f)).*sin(2*psi);
        Prr=Phh_mean(f).*(sin(psi)).^2 + Pdd_mean(f).*(cos(psi)).^2 - real(Phd_mean(f)).*sin(2*psi);
        
        [deltaSF_top(f,i),psi_max(f,i)]=max(Ptt./Prr);
        deltaSF(f,i)=real((deltaSF_top(f,i) - 1) / rms(tan(beta(f))));
    end
    %% Pz/Pg
    n_length=size(H_night,1);
    n_win=1024;
    n_ovrlp=0;
    n_fft=n_win;
    fs=1;
    
    % Exclude excessive points
    Z_unfilt0=Z_night(:,i); G_unfilt0=G_night(:,i); 
    
    excl_pts=rem(n_length,n_win);
    Z_unfilt=Z_unfilt0(1:end-excl_pts);  G_unfilt=G_unfilt0(1:end-excl_pts);
    
    % Power spectral density
    [Pzz,fff]=pwelch(Z_unfilt,hamming(n_win),n_ovrlp,n_fft,fs);
    [Pgg,fff]=pwelch(G_unfilt,hamming(n_win),n_ovrlp,n_fft,fs);
    
    % ULF extraction
    Pzz1=Pzz(bin_floor:bin_ceil,:);
    Pgg1=Pgg(bin_floor:bin_ceil,:);
     
    % ULF averaging
    for a=1:numel(f_rangeint)-1
        Pzz_band(:,a)=Pzz1(f_rangeidx(a):f_rangeidx(a+1)-1);
        Pgg_band(:,a)=Pgg1(f_rangeidx(a):f_rangeidx(a+1)-1);
        
        fff_band(:,a)=fff1(f_rangeidx(a):f_rangeidx(a+1)-1);
    end
    
    % Polarization ratio calculation
    ZG(:,i)=mean(abs(Pzz_band)./abs(Pgg_band),1);
    %% DepG
    Pgg_mean=mean(Pgg_band,1)';
    DepG0(:,i)=1./Pgg_mean.^2;
    if i>5
        DepG_back=mean(DepG0(:,i-5:i-1),2,'omitnan');
        DepG(:,i)=(DepG0(:,i)-DepG_back)./DepG_back;
    else
        DepG(:,i)=DepG0(:,i);
    end
    %% Direction estimation
    
    % Filtering
    H_pass1=filter1('bp',H_night(:,i),'fc',[f_rangeint(1) f_rangeint(end)],'fs',fs,'order',4);
    D_pass1=filter1('bp',D_night(:,i),'fc',[f_rangeint(1) f_rangeint(end)],'fs',fs,'order',4);
    
    % Azimuthal direction calculation
    H_hilb=hilbert(H_pass1); D_hilb=hilbert(D_pass1);
    
    PE_top=2*abs(H_hilb).*abs(D_hilb).*cos(angle(H_hilb)-angle(D_hilb));
    PE_bot=(abs(D_hilb)).^2-(abs(H_hilb)).^2;
    
    PE_theta=(atan2(PE_top,PE_bot))/2;
    PE_alpha1=pi-PE_theta;
    
    PE_alpha2=NaN(1,length(PE_alpha1));
    for b=1:numel(PE_alpha1)
        if PE_alpha1(b)==pi/2
            PE_alpha2(b)=3*pi/2;
        elseif PE_alpha1(b)>pi/2 && PE_alpha1(b)<pi
            PE_alpha2(b)=PE_alpha1(b)+pi;
        elseif PE_alpha1(b)==pi
            PE_alpha2(b)=0;
        elseif PE_alpha1(b)>pi && PE_alpha1(b)<3*pi/2
            PE_alpha2(b)=PE_alpha1(b)-pi;
        elseif PE_alpha1(b)==3*pi/2
            PE_alpha2(b)=pi/2;
        end
    end
    
    sum_HD=sqrt(H_pass0.^2+D_pass0.^2);
    mean_HD=sqrt(mean(H_pass.^2+D_pass.^2));
    
    for c=1:length(PE_alpha1)
        if sum_HD(c)>k_alpha*mean_HD
            alpha1(c,i)=PE_alpha1(c);
            alpha2(c,i)=PE_alpha2(c);
        end
    end
    
end
DepG(:,1:5)=NaN;
alpha3=vertcat(alpha1,alpha2);

%% Polar plotting parameters
EQ_plrscat=NaN(size(EQ_sel,1),5,days_num);
for j=1:days_num
    EQ_plrscat(:,1,j)=deg2rad(EQ_sel(:,8));                     % Angle in radian
    now_minus_EQ=datenum_vec(j)-round(EQ_sel(:,1));
    for k=1:size(EQ_sel,1)
        if now_minus_EQ(k)>=-9 && now_minus_EQ(k)<=0
            EQ_plrscat(k,2,j)=(10-abs(now_minus_EQ(k)))/10;     % Datenum (transpearancy)
        else
            EQ_plrscat(k,2,j)=0;
        end
        if now_minus_EQ(k)<=0 && now_minus_EQ(k)>=-29 
            EQ_plrscat(k,5,j)=0.5;                             % Datenum 2 (transperancy 2)
        else
            EQ_plrscat(k,5,j)=0;
        end
    end
    EQ_plrscat(:,3,j)=2.5.^EQ_sel(:,3);                         % Magnitude (circle size)
    EQ_plrscat(:,4,j)=100*(EQ_sel(:,2)/dis_max);                % Distance (distance from center)
end
%% Precursor parameters plotting

f1=figure(1);

stitle=sprintf('%s %s',stn,year);
suptitle(stitle)

sp1=subplot(4,1,1);
hold on
plot(hourly_vec,Dst,'b')
plot(hourly_vec,ones(length(hourly_vec),1)*(-50),'b--');
plot(hourly_vec,ap,'r')
plot(hourly_vec,ones(length(hourly_vec),1)*(50),'r--');
ylabel('Dst      ap')
hold off
yyaxis right
set(gca, 'YColor', 'k');
plot(datetime(EQ_sel(:,1),'ConvertFrom','datenum'),EQ_sel(:,4),'k^','MarkerFaceColor','m');
ylabel('K_{LS}')
xlim([min(days_vec) max(days_vec)]);
title('Geomagnetic (left) & local seismicity (right) indices')

sp2=subplot(4,1,2);
surf(days_vec,f_rangeint(1:end-1),deltaSF,'LineStyle','none','FaceColor','interp')
view([0,0,90])
xlim([min(days_vec) max(days_vec)]);
ylim([f_rangeint(1) f_rangeint(end-1)]);
ylabel('f (Hz)')
title('\DeltaS(f)')
colorbar('east')

sp3=subplot(4,1,3);
surf(days_vec,f_rangeint(1:end-1),ZG,'LineStyle','none','FaceColor','interp')
view([0,0,90])
xlim([min(days_vec) max(days_vec)]);
ylim([f_rangeint(1) f_rangeint(end-1)]);
ylabel('f (Hz)')
title('Z/G')
colorbar('east')

sp4=subplot(4,1,4);
surf(days_vec,f_rangeint(1:end-1),real(DepG),'LineStyle','none','FaceColor','interp')
view([0,0,90])
xlim([min(days_vec) max(days_vec)]);
ylim([f_rangeint(1) f_rangeint(end-1)]);
ylabel('f (Hz)')
title('DepG')
colorbar('east')

linkaxes([sp1,sp2,sp3,sp4],'x')
% %% Azimuthal beam polar plotting
% 
% f2=figure(2);
% 
% binedge=linspace(0,2*pi,37);
% polarscatter(NaN,NaN)
% set(gca,'ThetaDir','clockwise','ThetaZeroLocation','top');
% set(gcf, 'Position', get(0, 'Screensize'));
% hold on
% for i=1:days_num
% 
%     if any(~isnan(alpha3(:,i)))
%         histo=polarhistogram(alpha1(:,i),'BinEdges',binedge,'FaceColor','b');
%         hold on
%         polarhistogram(alpha2(:,i),'BinEdges',binedge,'FaceColor','b');
%         hold off
%         rlim([0 100]);
%         set(gca,'ThetaDir','clockwise','ThetaZeroLocation','top');
%         azim_exist=1;
%         drawnow;
%     else
%         polarhistogram(randn(1000,1),'BinEdges',binedge,'FaceAlpha',0,'EdgeAlpha',0);
%         rlim([0 100]);
%         set(gca,'ThetaDir','clockwise','ThetaZeroLocation','top');
%         azim_exist=0;
%         drawnow;
%     end
%     hold on
%     
%     if azim_exist==1
%        [azim_max(1,i),azim_ind]=max(histo.Values);
%        azim_max(2,i)=azim_max(1,i);
%        azim_maxangle(1,i)=mean([histo.BinEdges(azim_ind),histo.BinEdges(azim_ind+1)]);
%        if azim_maxangle(1,i)>pi/2 && azim_maxangle(1,i)<pi
%            azim_maxangle(2,i)=azim_maxangle(1,i)+pi; 
%        elseif azim_maxangle(1,i)>pi && azim_maxangle(1,i)<3*pi/2
%            azim_maxangle(2,i)=azim_maxangle(1,i)-pi;
%        end
%      
%     elseif azim_exist==0
%         azim_max(1,i)=NaN; azim_max(2,i)=NaN;
%         azim_maxangle(1,i)=NaN; azim_maxangle(2,i)=NaN;
%     end
%     
%     azim_title=sprintf('Date: %s',days_vec(i));
%     title(azim_title);
%     drawnow;
%     
%     for j=1:size(EQ_plrscat,1)
%         polarscatter(EQ_plrscat(j,1,i),EQ_plrscat(j,4,i),EQ_plrscat(j,3,i),'o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerFaceAlpha',EQ_plrscat(j,2,i));
%         set(gca,'ThetaDir','clockwise','ThetaZeroLocation','top');
%        drawnow;
%     end
%     rlim([0 100]);
%     hold off
%     pause(0.3)
% end
% 
% close(f2)
% 
% %% Single date azimuthal beam polar plotting
% azim_plot=input('Enter a date in [yyyy,mm,dd] format ');
% 
% f2=figure(2);
% 
% day_sel=find(datenum(days_vec)==datenum(azim_plot));
% 
% if any(~isnan(alpha1(:,day_sel)))
%     polarscatter(1,200,3,'o','MarkerFaceColor','r','MarkerEdgeColor','none');
%     hold on
%     legend('Upcoming EQs within 30 days','AutoUpdate','off','Location','southoutside')
%     polarhistogram(alpha1(:,day_sel),'BinEdges',binedge,'FaceColor','b');
%     polarhistogram(alpha2(:,day_sel),'BinEdges',binedge,'FaceColor','b');
%     rlim([0 100]);
%     for i=1:size(EQ_plrscat,1)
%         polarscatter(EQ_plrscat(i,1,day_sel),EQ_plrscat(i,4,day_sel),EQ_plrscat(i,3,day_sel),'o',...
%             'MarkerFaceColor','r','MarkerEdgeColor','none',...
%             'MarkerFaceAlpha',EQ_plrscat(i,5,day_sel));
%     end
%     hold off
%     set(gca,'ThetaDir','clockwise','ThetaZeroLocation','top');
% end
% title(sprintf('Azimuthal beam on %s',days_vec(day_sel)));
% %% Temporal evolution of maximum beam distribution
% 
% f3=figure(3);
% 
% title('Temporal evolution of maximum azimuthal distribution');
% plot(days_vec,rad2deg(azim_maxangle(1,:)),'Color',[1,0,0,0.3],'HandleVisibility','off');
% hold on
% k=1;
% for i=1:days_num
%     if ~isnan(azim_max(1,i))
%         sizemark=(100/max(azim_max(1,:)))*azim_max(1,i);   
%         scatter(days_vec(i),rad2deg(azim_maxangle(1,i)),azim_max(1,i),sizemark,'r','HandleVisibility','off');
%     end
% end
% for k=1:size(EQ_sel,1)
%     plot(datetime(EQ_sel(k,1),'ConvertFrom','datenum'),EQ_sel(k,8),'m^',...
%         'MarkerSize',EQ_sel(k,4)*(20/max(EQ_sel(:,4)))+5,'HandleVisibility','off');
% end
% plot(days_vec,rad2deg(azim_maxangle(2,:)),'Color',[0,0,1,0.3],'HandleVisibility','off');
% 
% for i=1:days_num
%     if ~isnan(azim_max(2,i))
%         sizemark=(100/max(azim_max(2,:)))*azim_max(2,i);   
%         scatter(days_vec(i),rad2deg(azim_maxangle(2,i)),azim_max(2,i),sizemark,'b','HandleVisibility','off');
%     end
% end
% for k=1:size(EQ_sel,1)
%     plot(datetime(EQ_sel(k,1),'ConvertFrom','datenum'),EQ_sel(k,8),'m^',...
%         'MarkerSize',EQ_sel(k,4)*(20/max(EQ_sel(:,4)))+5,'HandleVisibility','off' );
% end
% plot(days_vec(1),NaN,'ro',days_vec(1),NaN,'bo',days_vec(1),NaN,'m^');
% hold off
% title('Temporal evolution of maximum azimuthal beam');
% legend('\alpha_{1}','\alpha_{2}','\chi_{EQ}')
% ylim([0 360]); xlabel('Date'); ylabel(['Angle (',char(176),')']); yticks([0 90 180 270 360]);
% 
% 
%% Raw data plotting


if input('Plot original undenoised data? 1=Yes, 0=No ') plot_ori=1; else plot_ori=0; end

f4=figure(4);
title('Denoised and original data')
sec_vec=datetime(date_start1:1/86400:date_end1,'ConvertFrom','datenum','Format','dd/MM/yyyy HH:mm:ss');

sp1=subplot(4,1,1);
if plot_ori==1 plot(sec_vec,H_ori); end
hold on
plot(sec_vec,H_dn)
hold off
if plot_ori==1 legend('Original','Denoised'); else legend('Denoised'); end
ylabel('H (nT)')

sp2=subplot(4,1,2);
if plot_ori==1 plot(sec_vec,D_ori); end
hold on
plot(sec_vec,D_dn)
hold off
if plot_ori==1 legend('Original','Denoised'); else legend('Denoised'); end
ylabel('D (nT)')

sp3=subplot(4,1,3);
if plot_ori==1 plot(sec_vec,Z_ori); end
hold on
plot(sec_vec,Z_dn)
hold off
if plot_ori==1 legend('Original','Denoised'); else legend('Denoised'); end
ylabel('Z (nT)')

sp4=subplot(4,1,4);
plot(hourly_vec,ap)
hold on
plot(hourly_vec,Dst)
ylabel('ap   Dst')
xlim([min(hourly_vec) max(hourly_vec)])
hold off
linkaxes([sp1,sp2,sp3,sp4],'x')

% %% Single date raw data plotting
% 
% raw_plot=input('Enter a date in [yyyy,mm,dd] format ');
% 
% f5=figure(5);
% 
% day_sel1=find(datenum(days_vec)==datenum(raw_plot));
% time_sel=datetime(UT1m_night(:,day_sel1),'Format','dd/MM/yyyy HH:mm:ss','ConvertFrom','datenum');
% H_pass=filter1('bp',H_night(:,day_sel1),'fc',[f_rangeint(1) f_rangeint(end)],'fs',fs,'order',4);
% D_pass=filter1('bp',D_night(:,day_sel1),'fc',[f_rangeint(1) f_rangeint(end)],'fs',fs,'order',4);
% Z_pass=filter1('bp',Z_night(:,day_sel1),'fc',[f_rangeint(1) f_rangeint(end)],'fs',fs,'order',4);
% 
% sp1=subplot(3,1,1);
% plot(time_sel,H_night(:,day_sel1));
% ylabel('H_{raw}');
% 
% yyaxis right
% plot(time_sel,H_pass);
% ylabel('H_{ULF}');
% 
% sp2=subplot(3,1,2);
% plot(time_sel,D_night(:,day_sel1));
% ylabel('D_{raw}');
% yyaxis right
% plot(time_sel,D_pass);
% ylabel('D_{ULF}');
% 
% sp3=subplot(3,1,3);
% plot(time_sel,Z_night(:,day_sel1));
% ylabel('Z_{raw}');
% yyaxis right
% plot(time_sel,Z_pass);
% ylabel('Z_{ULF}');
% 
% linkaxes([sp1,sp2,sp3],'x')
%% Saving figures

for i=1:5
    input_figname0=input('Name of figure you want to save? Insert number only. 0=Finish saving  ','s');
    input_figname=['f',input_figname0];
    
    if strcmp(input_figname0,'0') disp('Saving finished.'); break; end
    ddMMyy1=upper(string(datetime(min(days_vec),'Format','ddMMyy')));
    ddMMyy2=upper(string(datetime(max(days_vec),'Format','ddMMyy')));
    
    today=strcat(char(datetime('today','Format','dd-MM-yyyy')),'\Adib_103\');
    path1='D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\';
    mkdir(fullfile(path1,today));
    path2=strcat('D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\',today);
    
    figname1=char(strcat(stn,'_',ddMMyy1,'-',ddMMyy2));
    
    for i=1:100
        if exist(strcat(path2,figname1,sprintf('-%d',i),'.tif'))~=2
            figname1ext=strcat(figname1,char(sprintf('-%d',i)));
            saveas(eval(input_figname),fullfile(path2,figname1ext),'tiff');
            figpath=strcat(path2,figname1ext,'.tif');
            break;
        end
    end
    
end

