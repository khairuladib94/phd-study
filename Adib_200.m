%% Purpose
% Observe Pz/Pg in 9 frequency bands, e.g. 0.01 - 0.02, 0.02 - 0.03 .. 0.09
% - 0.10 hz
% Estimate angle using polarization ellipse method as used by Schekotov
% The frequency band to bandpass filter horizontal components to calculate
% azimuthal angle can be chosen, e.g. 0.01 - 0.02, 0.02 - 0.03 .. 0.09 - 0.10 hz
% 2007, 2008, 2015 and Ohta 2013
% Altered from Adib_104
%% Load initial variable
clc; close all; clearvars -except EQ_intrst;
load VARIABLES_WORLD
open EQ_intrst
%% Customization
EQ_num=43;          %Select earthquake of interest as listed in EQ_intrst #

mag_min=5.0;        %Minimum magnitude of earthquakes to be considered
dis_max=300;        %Maximum epicentral distance from the station
depth_max=200;

LT_sel=1;
LT_start=[22,23,00,01,02]; LT_start=LT_start(LT_sel);
LT_end=  [02,03,04,05,06]; LT_end=LT_end(LT_sel);

f_rangeint=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10];

k_Thres=7; N_medfil=1;  k_repeat=4;       %Noise removal parameters

k_alpha=3.5;                              %Minimal SNR coefficient
            
normalize_ZG=1;                           %1 to normalize ZG

%Frequency range to use for bandpass filter for azimuth direction.
%1=0.01-0.02,2=0.02-0.03...9=0.09-0.10 Hz. Can enter multiple values.
f_prec=3;

daysafter_prec=30;                  %Number of days to show upcoming earthquakes
remove_disturbed=0;                 %Keep (0) or remove (1) data during disturbed days
%% Place and date
stn0=EQ_intrst{EQ_num,1}; stn=stn0(1:3);
EQ_date0=EQ_intrst{EQ_num,3}; 
if numel(EQ_date0)<14 
    for i=14:-1:11 EQ_date0(i)=EQ_date0(i-1); end 
    EQ_date0(10)='0'; 
end
EQ_date1=[str2num(strcat('20',EQ_date0(7:8))),str2num(EQ_date0(4:5)),str2num(EQ_date0(1:2)),str2num(EQ_date0(10:11)),str2num(EQ_date0(13:14)),0];
EQ_date=datenum(EQ_date1(1:3));
date_start0=EQ_date-45;
date_end0=EQ_date+15;
date_start=datevec(date_start0);
date_end=datevec(date_end0);
%% Time period
date_start1=datenum(horzcat(date_start(1:3),[0,0,0]));
date_end1=datenum(horzcat(date_end(1:3),[23,59,59]));
datenum_start=datenum(date_start);
datenum_end=datenum(date_end);
days_num=datenum_end-datenum_start+1;
days_vec=datetime(datevec(datenum_start+0.5:1:datenum_end+0.5),'Format','dd/MM/yyyy HH');

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
EQ_sel=NaN(1,9);
for j=1:size(EQ_table,1)
    if  EQ_table(j,10)>=mag_min
        EQ_latlon=[EQ_table(j,7:8)];
        EQ_dis=deg2km(distance('gc',stn_latlon,EQ_latlon));
        EQ_time=datenum(EQ_table(j,1:6));
        EQ_mag=EQ_table(j,10);
        EQ_Ks=(10^(0.75*EQ_mag))/(EQ_dis+100);
        EQ_depth=EQ_table(j,9);
        EQ_angle=azimuth('gc',stn_latlon(1:2),EQ_latlon(1:2),'degree');
        if 0.025*EQ_dis<EQ_mag-4.5
            EQ_detect=1;
        else
            EQ_detect=0;
        end
        
        if EQ_dis<=dis_max && EQ_depth<=depth_max
            EQ_sel(i,1)=EQ_time;
            EQ_sel(i,2)=EQ_dis;
            EQ_sel(i,3)=EQ_mag;
            EQ_sel(i,4)=EQ_Ks;
            EQ_sel(i,5)=EQ_depth;
            EQ_sel(i,6:7)=EQ_latlon;
            EQ_sel(i,8)=EQ_angle;
            EQ_sel(i,9)=EQ_detect;
            i=i+1;
        end
    end
end

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
%% Remove data -1, on and +1 during disturbed days
ap_daily=reshape(ap,[],days_num);
Dst_daily=reshape(Dst,[],days_num);
if remove_disturbed==1
    for i=1:days_num
        if any(ap_daily(:,i)>50)||any(ap_daily(:,i)<-50)||any(Dst_daily(:,i)>50)||any(Dst_daily(:,i)<-50)
            H_night(:,i)=NaN;
            D_night(:,i)=NaN;
            Z_night(:,i)=NaN;
            G_night(:,i)=NaN;
        end
    end
end
%% Polar plotting parameters
EQ_plrscat=NaN(size(EQ_sel,1),6,days_num);
for j=1:days_num
    EQ_plrscat(:,1,j)=deg2rad(EQ_sel(:,8));                     % Angle in radian
    now_minus_EQ=datenum_vec(j)-round(EQ_sel(:,1));
    for k=1:size(EQ_sel,1)
        if now_minus_EQ(k)>=-9 && now_minus_EQ(k)<=0
            EQ_plrscat(k,2,j)=(10-abs(now_minus_EQ(k)))/10;     % Datenum (transparency)
        else
            EQ_plrscat(k,2,j)=0;
        end
        if now_minus_EQ(k)<=0 && now_minus_EQ(k)>=-(daysafter_prec-1)
            EQ_plrscat(k,5,j)=0.5;                             % Datenum 2 (transparency 2)
            EQ_plrscat(k,6,j)=EQ_sel(k,9);
        else
            EQ_plrscat(k,5,j)=0;
            EQ_plrscat(k,6,j)=0;
        end
    end
    EQ_plrscat(:,3,j)=2.5.^EQ_sel(:,3);                         % Magnitude (circle size)
    EQ_plrscat(:,4,j)=100*(EQ_sel(:,2)/dis_max);                % Distance (distance from center)
    
end
%% Precursor detection

alpha1=NaN(size(H_night,1),days_num);
alpha2=NaN(size(H_night,1),days_num);

for i=1:days_num
    
    % Skip if data contains nans
    if any(isnan(Z_night(:,i))) || any(isnan(G_night(:,i)))
        ZZ(:,i)=NaN(1,9);
        GG(:,i)=NaN(1,9);
        continue
    end
    
    % Power spectral density
    n_length=size(H_night,1);
    n_win=n_length;
    n_ovrlp=0;
    n_fft=n_win;
    fs=1;
    
    [Pzz,fff]=periodogram(Z_night(:,i),hamming(n_win),n_fft,fs);
    [Pgg,fff]=periodogram(G_night(:,i),hamming(n_win),n_fft,fs);
    
    bin_floor=find(round(fff,4)==f_rangeint(1)); bin_ceil=find(round(fff,4)==f_rangeint(end));
    
    % ULF extraction
    Pzz1=Pzz(bin_floor:bin_ceil,:);
    Pgg1=Pgg(bin_floor:bin_ceil,:);
    fff1=fff(bin_floor:bin_ceil,:);
    
    % ULF averaging
    for b=1:numel(f_rangeint) f_rangeidx(b)=find(round(fff1,4)==f_rangeint(b)); end
    
    for a=1:numel(f_rangeint)-1
        Pzz_band(:,a)=Pzz1(f_rangeidx(a):f_rangeidx(a+1)-1);
        Pgg_band(:,a)=Pgg1(f_rangeidx(a):f_rangeidx(a+1)-1);
        
        fff_band(:,a)=fff1(f_rangeidx(a):f_rangeidx(a+1)-1);
    end
    
    % Polarization ratio calculation
    ZZ(:,i)=real(mean(Pzz_band));
    GG(:,i)=real(mean(Pgg_band));
    
end 

if normalize_ZG==1
    ZZ_norm=normalize(ZZ,2,'range',[1 3]);
    GG_norm=normalize(GG,2,'range',[1 2]);
    ZG=ZZ_norm./GG_norm;                     % Normalized
else
    ZG=ZZ./GG;                               % Not normalized
end
%% Direction estimation & plotting

for j=1:numel(f_prec)
    alpha1=NaN(size(H_night,1),days_num);
    alpha2=NaN(size(H_night,1),days_num);
    
    for i=1:days_num
        
        if any(isnan(H_night(:,i))) || any(isnan(D_night(:,i)))
           continue; 
        end
        
        % Filtering
        H_pass=filter1('bp',H_night(:,i),'fc',[f_rangeint(f_prec(j)) f_rangeint(f_prec(j)+1)],'fs',fs,'order',4);
        D_pass=filter1('bp',D_night(:,i),'fc',[f_rangeint(f_prec(j)) f_rangeint(f_prec(j)+1)],'fs',fs,'order',4);
        
        % Azimuthal direction calculation
        H_hilb=hilbert(H_pass); D_hilb=hilbert(D_pass);
        
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
        
        sum_HD=sqrt(H_pass.^2+D_pass.^2);
        mean_HD=sqrt(mean(H_pass.^2+D_pass.^2));
        
        for c=1:length(PE_alpha1)
            if sum_HD(c)>k_alpha*mean_HD
                alpha1(c,i)=PE_alpha1(c);
                alpha2(c,i)=PE_alpha2(c);
            end
        end
        
    end
    
    alpha3=vertcat(alpha1,alpha2);
    %% Temporal evolution of maximum beam distribution
    
    f3=figure(3);
%     if j==1
%         stitle=sprintf('%s %s',stn,year);
%         sptitle=suptitle(stitle);
%         stitle_pos=get(sptitle,'position');
%         stitle_pos=[stitle_pos(1) stitle_pos(2)+0.02 stitle_pos(3)];
%         set(sptitle,'position',stitle_pos);
%     end
    
    subplot(numel(f_prec),1,j)
    
    binedge=linspace(0,2*pi,37);
    for i=1:days_num
        if any(~isnan(alpha3(:,i)))
            set(0,'DefaultFigureVisible','off');
            histo=polarhistogram(alpha1(:,i),'BinEdges',binedge,'FaceColor','b');
            [azim_max(1,i),azim_ind]=max(histo.Values);
            azim_max(2,i)=azim_max(1,i);
            azim_maxangle(1,i)=mean([histo.BinEdges(azim_ind),histo.BinEdges(azim_ind+1)]);
            if azim_maxangle(1,i)>pi/2 && azim_maxangle(1,i)<pi
                azim_maxangle(2,i)=azim_maxangle(1,i)+pi;
            elseif azim_maxangle(1,i)>pi && azim_maxangle(1,i)<3*pi/2
                azim_maxangle(2,i)=azim_maxangle(1,i)-pi;
            end
        else
            azim_max(1,i)=NaN; azim_max(2,i)=NaN;
            azim_maxangle(1,i)=NaN; azim_maxangle(2,i)=NaN;
        end
    end
    
    set(0,'DefaultFigureVisible','on')
    
    if j==1 title('Temporal evolution of maximum azimuthal distribution'); end
    plot(days_vec,rad2deg(azim_maxangle(1,:)),'Color',[1,0,0,0.3],'HandleVisibility','off');
    hold on
    % alpha 1
    k=1;
    for i=1:days_num
        if ~isnan(azim_max(1,i))
            sizemark=(100/max(azim_max(1,:)))*azim_max(1,i);
            scatter(days_vec(i),rad2deg(azim_maxangle(1,i)),azim_max(1,i),sizemark,'r','HandleVisibility','off');
        end
    end
    % alpha 2
    for i=1:days_num
        if ~isnan(azim_max(2,i))
            sizemark=(100/max(azim_max(2,:)))*azim_max(2,i);
            scatter(days_vec(i),rad2deg(azim_maxangle(2,i)),azim_max(2,i),sizemark,'b','HandleVisibility','off');
        end
    end
    plot(days_vec,rad2deg(azim_maxangle(2,:)),'Color',[0,0,1,0.3],'HandleVisibility','off');
    % EQ
    for k=1:size(EQ_sel,1)
        if EQ_sel(k,9)==1
            plot(datetime(EQ_sel(k,1),'ConvertFrom','datenum'),EQ_sel(k,8),'m^','LineWidth',2,...
                'MarkerSize',EQ_sel(k,4)*(20/max(EQ_sel(:,4)))+5,'HandleVisibility','off');
        else
            plot(datetime(EQ_sel(k,1),'ConvertFrom','datenum'),EQ_sel(k,8),'m^','MarkerFaceColor','none',...
                'MarkerSize',EQ_sel(k,4)*(20/max(EQ_sel(:,4)))+5,'HandleVisibility','off');
        end
    end
    line([datetime(EQ_date1) datetime(EQ_date1)],[-1e5 1e5],'Color','k','LineWidth',0.05,'HandleVisibility','off');
    plot(days_vec(1),NaN,'ro',days_vec(1),NaN,'bo',days_vec(1),NaN,'m^');
    hold off
    title(sprintf('Temporal evolution of maximum azimuthal beam (%cf_{BP}=%.2f-%.2f Hz)',char(916),f_rangeint(f_prec(j)),f_rangeint(f_prec(j)+1)));
    legend('\alpha_{1}','\alpha_{2}','\chi_{EQ}')
    ylim([0 360]); xlabel('Date'); ylabel(['Angle (',char(176),')']); yticks([0 90 180 270 360]);
end
%% Precursor parameters plotting
f1=figure(1);

% stitle=sprintf('%s %s',stn,year);
% sptitle=suptitle(stitle);
% stitle_pos=get(sptitle,'position');
% stitle_pos=[stitle_pos(1) stitle_pos(2)+0.02 stitle_pos(3)];
% set(sptitle,'position',stitle_pos);

colpal=lines(3);

sp1=subplot(4,1,1);
hold on
for k=1:days_num
    if any(ap_daily(:,k)>50) || any(Dst_daily(:,k)<-50)
        xpatch=[k-1 k k k-1];
        ypatch=[-1e+5 -1e+5 1e+5 1e+5];
        patch(xpatch,ypatch,'r','FaceAlpha',0.1,'EdgeAlpha',0,'HandleVisibility','off');
    end
end
plot(hourly_vec,Dst,'b')
plot(hourly_vec,ones(length(hourly_vec),1)*(-50),'b--');
plot(hourly_vec,ap,'r')
plot(hourly_vec,ones(length(hourly_vec),1)*(50),'r--');
    ylabel('Dst      ap'); ylim([min(Dst) max(ap)]);
hold off
yyaxis right
set(gca, 'YColor', 'k');
hold on
plot(datetime(EQ_sel(find(EQ_sel(:,9)==1),1),'ConvertFrom','datenum'),EQ_sel(find(EQ_sel(:,9)==1),4),'m^','MarkerFaceColor','m');
plot(datetime(EQ_sel(find(EQ_sel(:,9)==0),1),'ConvertFrom','datenum'),EQ_sel(find(EQ_sel(:,9)==0),4),'m^','MarkerFaceColor','none');
for i=1:size(EQ_sel,1)
    if EQ_sel(i,4)>=min(maxk(EQ_sel(:,4),5)) || EQ_sel(i,9)==1
        label=sprintf('  M%0.1f, h=%.0fKM, d=%.0fKM',EQ_sel(i,3),EQ_sel(i,5),EQ_sel(i,2));
        text(datetime(EQ_sel(i,1),'ConvertFrom','datenum'),EQ_sel(i,4),label,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',7.5);
    end
end
ylim([0 max(EQ_sel(:,4))]);
line([datetime(EQ_date1) datetime(EQ_date1)],[-1e5 1e5],'Color','k');
hold off

ylabel('K_{LS}')
xlim([min(days_vec) max(days_vec)]);
title('ap & Dst (left) and K_{LS} (right)')

sp2=subplot(4,1,2);
hold on
for k=1:days_num
    if any(ap_daily(:,k)>50) || any(Dst_daily(:,k)<-50)
        xpatch=[k-1 k k k-1];
        ypatch=[-1e+5 -1e+5 1e+5 1e+5];
        patch(xpatch,ypatch,'r','FaceAlpha',0.1,'EdgeAlpha',0,'HandleVisibility','off');
    end
end
for i=1:3
    plot(days_vec,ZG(i,:),'-o','Color',colpal(i,:));
    ZG_thres=nanmean(ZG(i,:))+2*nanstd(ZG(i,:));
    plot(days_vec,ones(1,numel(days_vec))*ZG_thres,'--','Color',colpal(i,:),'HandleVisibility','off')
end
line([datetime(EQ_date1) datetime(EQ_date1)],[-1e5 1e5],'Color','k');
hold off
xlim([min(days_vec) max(days_vec)]); ylim([0 3.0]);
ylabel('a.u.')
legend('0.01-0.02','0.02-0.03','0.03-0.04');
title(sprintf('Z/G during %cLT=%02d-%02d (LT_{sel}=%d)',char(916),LT_start,LT_end,LT_sel))

sp3=subplot(4,1,3);
hold on
for k=1:days_num
    if any(ap_daily(:,k)>50) || any(Dst_daily(:,k)<-50)
        xpatch=[k-1 k k k-1];
        ypatch=[-1e+5 -1e+5 1e+5 1e+5];
        patch(xpatch,ypatch,'r','FaceAlpha',0.1,'EdgeAlpha',0,'HandleVisibility','off');
    end
end
for i=4:6
    plot(days_vec,ZG(i,:),'-o','Color',colpal(i-3,:));
    ZG_thres=nanmean(ZG(i,:))+2*nanstd(ZG(i,:));
    plot(days_vec,ones(1,numel(days_vec))*ZG_thres,'--','Color',colpal(i-3,:),'HandleVisibility','off')
end
line([datetime(EQ_date1) datetime(EQ_date1)],[-1e5 1e5],'Color','k');
hold off
xlim([min(days_vec) max(days_vec)]); ylim([0 3.0]);
ylabel('a.u.')
legend('0.04-0.05','0.05-0.06','0.06-0.07');

sp4=subplot(4,1,4);
hold on
for k=1:days_num
    if any(ap_daily(:,k)>50) || any(Dst_daily(:,k)<-50)
        xpatch=[k-1 k k k-1];
        ypatch=[-1e+5 -1e+5 1e+5 1e+5];
        patch(xpatch,ypatch,'r','FaceAlpha',0.1,'EdgeAlpha',0,'HandleVisibility','off');
    end
end
for i=7:9
    plot(days_vec,ZG(i,:),'-o','Color',colpal(i-6,:));
    ZG_thres=nanmean(ZG(i,:))+2*nanstd(ZG(i,:));
    plot(days_vec,ones(1,numel(days_vec))*ZG_thres,'--','Color',colpal(i-6,:),'HandleVisibility','off')
end
line([datetime(EQ_date1) datetime(EQ_date1)],[-1e5 1e5],'Color','k');
hold off
xlim([min(days_vec) max(days_vec)]); ylim([0 3.0]);
ylabel('a.u.')
legend('0.07-0.08','0.08-0.09','0.09-0.10');

linkaxes([sp1,sp2,sp3,sp4],'x')
linkaxes([sp2,sp3,sp4],'y')
%% Stop script running
return
%% Single date azimuthal beam polar plotting
azim_plot=input('Enter a date in [yyyy,mm,dd] format ');
day_sel=find(floor(datenum(days_vec))==datenum(azim_plot));

f2=figure(2);

set(0,'DefaultFigureVisible','off');
histo1=polarhistogram(alpha3(:,day_sel),'BinEdges',binedge,'FaceColor','b');
azim_total=histo1.Values;

set(0,'DefaultFigureVisible','on');
if any(~isnan(alpha1(:,day_sel)))
    polarscatter(1,200,3,'o','MarkerFaceColor','r','MarkerEdgeColor','none');
    legend(sprintf('Upcoming EQs within %d days',daysafter_prec),'AutoUpdate','off','Location','southoutside');
    hold on
    for j=1:numel(azim_total)
        polarhistogram(mean(binedge(j:j+1))*ones(100,1),'BinEdges',binedge,...
            'FaceColor','b','FaceAlpha',(azim_total(j)-min(azim_total))/range(azim_total),'EdgeColor','none');
    end
    rlim([0 100]);
    for i=1:size(EQ_plrscat,1)
        polarscatter(EQ_plrscat(i,1,day_sel),EQ_plrscat(i,4,day_sel),EQ_plrscat(i,3,day_sel),'o',...
            'MarkerFaceColor','r','MarkerEdgeColor','r',...
            'MarkerFaceAlpha',EQ_plrscat(i,5,day_sel),'MarkerEdgeAlpha',EQ_plrscat(i,6,day_sel));
    end
    hold off
    set(gca,'ThetaDir','clockwise','ThetaZeroLocation','top');
end
title(sprintf('Azimuthal beam on %s (%cf_{BP}=%.2f-%.2f Hz)',days_vec(day_sel),char(916),f_rangeint(f_prec),f_rangeint(f_prec+1)));
%% Raw data plotting

stitle=sprintf('%s %s',stn,year);
suptitle(stitle);

f4=figure(4);

stitle=sprintf('%s %s',stn,year);
sptitle=suptitle(stitle);
stitle_pos=get(sptitle,'position');
stitle_pos=[stitle_pos(1) stitle_pos(2)+0.02 stitle_pos(3)];
set(sptitle,'position',stitle_pos);

title('Denoised and original data')
sec_vec=datetime(date_start1:1/86400:date_end1,'ConvertFrom','datenum','Format','dd/MM/yyyy HH:mm:ss');

sp1=subplot(5,1,1);
plot(sec_vec,H_ori);
hold on
plot(sec_vec,H_dn)
hold off
legend('Original','Denoised');
ylabel('H (nT)')
xlim([sec_vec(1) sec_vec(end)]); ylim([min(H_dn) max(H_dn)]);
title(sprintf('Raw data'));

sp2=subplot(5,1,2);
plot(sec_vec,D_ori);
hold on
plot(sec_vec,D_dn)
hold off
legend('Original','Denoised');
ylabel('D (nT)')
xlim([sec_vec(1) sec_vec(end)]); ylim([min(D_dn) max(D_dn)]);

sp3=subplot(5,1,3);
plot(sec_vec,Z_ori);
hold on
plot(sec_vec,Z_dn)
hold off
legend('Original','Denoised');
ylabel('Z (nT)')
xlim([sec_vec(1) sec_vec(end)]); ylim([min(Z_dn) max(Z_dn)]);

sp4=subplot(5,1,4);
plot(sec_vec,G_ori);
hold on
plot(sec_vec,G_dn)
hold off
legend('Original','Denoised');
ylabel('G (nT)')
xlim([sec_vec(1) sec_vec(end)]); ylim([min(G_dn) max(G_dn)]);

sp5=subplot(5,1,5);
plot(hourly_vec,Dst,'b')
hold on
plot(hourly_vec,ap,'r')
plot(hourly_vec,ones(length(hourly_vec),1)*(-50),'b--');
plot(hourly_vec,ones(length(hourly_vec),1)*(50),'r--');
ylabel('Dst   ap')
xlim([min(hourly_vec) max(hourly_vec)])
hold off
linkaxes([sp1,sp2,sp3,sp4,sp5],'x')

%% Single date raw data plotting

raw_plot=input('Enter a date in [yyyy,mm,dd] format ');

f5=figure(5);

% stitle=sprintf('%s %s',stn,year);
% sptitle=suptitle(stitle);
% stitle_pos=get(sptitle,'position');
% stitle_pos=[stitle_pos(1) stitle_pos(2)+0.02 stitle_pos(3)];
% set(sptitle,'position',stitle_pos);

day_sel1=find(floor(datenum(days_vec))==datenum(raw_plot));
time_sel=datetime(UT1m_night(:,day_sel1),'Format','dd/MM/yyyy HH:mm:ss','ConvertFrom','datenum');
Z_pass=filter1('bp',Z_night(:,day_sel1),'fc',[f_rangeint(f_prec) f_rangeint(f_prec+1)],'fs',fs,'order',4);
G_pass=filter1('bp',G_night(:,day_sel1),'fc',[f_rangeint(f_prec) f_rangeint(f_prec+1)],'fs',fs,'order',4);

ap_secly0=ap_daily(t_start:t_start+LT_length-1,day_sel1);
Dst_secly0=Dst_daily(t_start:t_start+LT_length-1,day_sel1);
k=1;
for i=0:3600:10800
    ap_secly(i+1:i+3600)=ap_secly0(k);
    Dst_secly(i+1:i+3600)=Dst_secly0(k);
    k=k+1;
end

sp_1=subplot(3,1,1);
plot(time_sel,Z_night(:,day_sel1));
ylabel('Z_{raw}');
yyaxis right
plot(time_sel,Z_pass);
ylabel('Z_{ULF}');

sp_2=subplot(3,1,2);
plot(time_sel,G_night(:,day_sel1));
ylabel('G_{raw}');
yyaxis right
plot(time_sel,G_pass);
ylabel('G_{ULF}');

sp_3=subplot(3,1,3);
plot(time_sel,ap_secly,'r',time_sel,Dst_secly,'b');
hold on
plot(time_sel,50*ones(1,numel(time_sel)),'r--',time_sel,-50*ones(1,numel(time_sel)),'b--');
hold off
ylabel('Dst           ap')
xlabel('Time')

linkaxes([sp_1,sp_2,sp_3],'x')

f6=figure(6);
% stitle=sprintf('%s %s',stn,year);
% sptitle=suptitle(stitle);
% stitle_pos=get(sptitle,'position');
% stitle_pos=[stitle_pos(1) stitle_pos(2)+0.02 stitle_pos(3)];
% set(sptitle,'position',stitle_pos);

[Pzz,fff]=periodogram(Z_night(:,day_sel1),hamming(n_win),n_fft,fs);
[Pgg,fff]=periodogram(G_night(:,day_sel1),hamming(n_win),n_fft,fs);

sp_4=subplot(2,1,1);
plot(fff,Pzz);
xlim([0.01 0.1]);
title(sprintf('Periodograms of Z and G on %s',days_vec(day_sel1)));
xlabel('f (Hz)'); ylabel('Power of Z');

sp_5=subplot(2,1,2);
plot(fff,Pgg);
xlim([0.01 0.1]); 
xlabel('f (Hz)'); ylabel('Power of G');

linkaxes([sp_4,sp_5],'xy');

%% ZZ, GG and ZG plotting over frequency
raw_plot=input('Enter a date in [yyyy,mm,dd] format ');
day_sel1=find(floor(datenum(days_vec))==datenum(raw_plot));

f7=figure(7);

% stitle=sprintf('%s %s',stn,year);
% sptitle=suptitle(stitle);
% stitle_pos=get(sptitle,'position');
% stitle_pos=[stitle_pos(1) stitle_pos(2)+0.02 stitle_pos(3)];
% set(sptitle,'position',stitle_pos);[]

bar(f_rangeint(1:end-1)+0.005,ZZ_norm(:,day_sel1),'FaceAlpha',0.5)
hold on
bar(f_rangeint(1:end-1)+0.005,GG_norm(:,day_sel1),'FaceAlpha',0.5)
plot(f_rangeint(1:end-1)+0.005,ZG(:,day_sel1),'-o')
hold off
xlim([0.01 0.1]); ylim([0 3]);
title(sprintf('Normalized Z, G and Z/G on %s',days_vec(day_sel1)));
legend('ZZ_{norm}','GG_{norm}','Z/G')
%% Saving figures
ddMMyy1=upper(string(datetime(min(days_vec),'Format','ddMMyy')));
ddMMyy2=upper(string(datetime(max(days_vec),'Format','ddMMyy')));

today=strcat(char(datetime('today','Format','dd-MM-yyyy')),'\Adib_200\');
path1='D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\';
mkdir(fullfile(path1,today));
path2=strcat('D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\',today);

figname1=char(strcat(stn,'_',ddMMyy1,'-',ddMMyy2));

input_figname0=input('Name of figure you want to save? Insert number in [ ] or 99=All  ');
if input_figname0==99
    clear input_figname0
    all_fig0=findobj('Type','figure');
    for j=1:numel(all_fig0)
        input_figname0(j)=all_fig0(j).Number;
    end
end

for l=1:numel(input_figname0)
    input_figname=strcat('f',string(input_figname0(l)));
    set(eval(input_figname), 'Position', get(0, 'Screensize'));
    for k=1:100
        if exist(strcat(path2,figname1,sprintf('-%d',k),'.tif'))~=2
            figname1ext=strcat(figname1,char(sprintf('-%d',k)));
            saveas(eval(input_figname),fullfile(path2,figname1ext),'tiff');
            figpath=strcat(path2,figname1ext,'.tif');
            break;
        end
    end
    close(eval(input_figname))
end
disp('Saving completed');
%% Azimuthal beam polar plotting
if input('Show azimuthal beam plotting? 1=Yes, 0=No')
    f2=figure(2);
    
    binedge=linspace(0,2*pi,37);
    polarscatter(NaN,NaN)
    set(gca,'ThetaDir','clockwise','ThetaZeroLocation','top');
    set(gcf, 'Position', get(0, 'Screensize'));
    hold on
    for i=1:days_num
        
        if any(~isnan(alpha3(:,i)))
            histo=polarhistogram(alpha1(:,i),'BinEdges',binedge,'FaceColor','b');
            hold on
            polarhistogram(alpha2(:,i),'BinEdges',binedge,'FaceColor','b');
            hold off
            rlim([0 100]);
            set(gca,'ThetaDir','clockwise','ThetaZeroLocation','top');
            azim_exist=1;
            drawnow;
        else
            polarhistogram(randn(1000,1),'BinEdges',binedge,'FaceAlpha',0,'EdgeAlpha',0);
            rlim([0 100]);
            set(gca,'ThetaDir','clockwise','ThetaZeroLocation','top');
            azim_exist=0;
            drawnow;
        end
        hold on
        
        if azim_exist==1
            [azim_max(1,i),azim_ind]=max(histo.Values);
            azim_max(2,i)=azim_max(1,i);
            azim_maxangle(1,i)=mean([histo.BinEdges(azim_ind),histo.BinEdges(azim_ind+1)]);
            if azim_maxangle(1,i)>pi/2 && azim_maxangle(1,i)<pi
                azim_maxangle(2,i)=azim_maxangle(1,i)+pi;
            elseif azim_maxangle(1,i)>pi && azim_maxangle(1,i)<3*pi/2
                azim_maxangle(2,i)=azim_maxangle(1,i)-pi;
            end
            
        elseif azim_exist==0
            azim_max(1,i)=NaN; azim_max(2,i)=NaN;
            azim_maxangle(1,i)=NaN; azim_maxangle(2,i)=NaN;
        end
        
        azim_title=sprintf('Date: %s',days_vec(i));
        title(azim_title);
        drawnow;
        
        for j=1:size(EQ_plrscat,1)
            polarscatter(EQ_plrscat(j,1,i),EQ_plrscat(j,4,i),EQ_plrscat(j,3,i),'o','MarkerFaceColor','r','MarkerEdgeColor','none','MarkerFaceAlpha',EQ_plrscat(j,2,i));
            set(gca,'ThetaDir','clockwise','ThetaZeroLocation','top');
            drawnow;
        end
        rlim([0 100]);
        hold off
        pause(0.3)
    end
    
    close(f2)
end
%% For thesis - precursor (single frequency)
f10=figure(10);

% stitle=sprintf('       %s %s',stn,year);
% sptitle=suptitle(stitle);
% stitle_pos=get(sptitle,'position');
% stitle_pos=[stitle_pos(1) stitle_pos(2)+0.02 stitle_pos(3)];
% set(sptitle,'position',stitle_pos);

sp1=subplot(2,1,1);

date_EQ1=datevec(EQ_sel(:,1));
date_EQ2=date_EQ1;
date_EQ2(:,6)=0;
EQ_ID=find(datenum(EQ_date1)==datenum(date_EQ2));

hold on
for k=1:days_num
    if any(ap_daily(:,k)>50) || any(Dst_daily(:,k)<-50)
        xpatch=[k-1 k k k-1];
        ypatch=[-1e+5 -1e+5 1e+5 1e+5];
        patch(xpatch,ypatch,'r','FaceAlpha',0.1,'EdgeAlpha',0,'HandleVisibility','off');
    end
end
plot(hourly_vec,Dst,'b')
plot(hourly_vec,ones(length(hourly_vec),1)*(-50),'b--','HandleVisibility','off');
plot(hourly_vec,ap,'r')
plot(hourly_vec,ones(length(hourly_vec),1)*(50),'r--','HandleVisibility','off');
    ylabel('ap & Dst (nT)'); ylim([min(Dst) max(ap)]);
hold off
legend('Dst','ap','Location','northwest')
yyaxis right
set(gca, 'YColor', 'k');
hold on
for i=1:size(EQ_sel,1)
    if i==EQ_ID
        plot(datetime(EQ_sel(i,1),'ConvertFrom','datenum'),EQ_sel(i,4),'m^','MarkerFaceColor','m','HandleVisibility','off');
        label=sprintf('  %s LT \n  M%0.1f, h=%.0fKM, d=%.0fKM',datetime(EQ_date1),EQ_sel(i,3),EQ_sel(i,5),EQ_sel(i,2));
        text(datetime(EQ_sel(i,1),'ConvertFrom','datenum'),EQ_sel(i,4),label,'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',9);
    else
        plot(datetime(EQ_sel(i,1),'ConvertFrom','datenum'),EQ_sel(i,4),'m^','MarkerFaceColor','none','HandleVisibility','off');
    end
end
ylim([0 max(EQ_sel(:,4))]);
line([datetime(EQ_date1) datetime(EQ_date1)],[-1e5 1e5],'Color','k','HandleVisibility','off');
hold off
ylabel('K_{LS}'); 
xlim([min(days_vec) max(days_vec)]);
title('ap & Dst (left) and K_{LS} (right)');
xticks([]);

sp2=subplot(2,1,2);
hold on
for k=1:days_num
    if any(ap_daily(:,k)>50) || any(Dst_daily(:,k)<-50)
        xpatch=[k-1 k k k-1];
        ypatch=[-1e+5 -1e+5 1e+5 1e+5];
        patch(xpatch,ypatch,'r','FaceAlpha',0.1,'EdgeAlpha',0,'HandleVisibility','off');
    end
end
plot(days_vec,ZG(f_prec,:),'-o','Color','k');
ZG_thres=nanmean(ZG(f_prec,:))+2*nanstd(ZG(f_prec,:));
plot(days_vec,ones(1,numel(days_vec))*ZG_thres,'--','Color','k','HandleVisibility','off')
line([datetime(EQ_date1) datetime(EQ_date1)],[-1e5 1e5],'Color','k');
hold off
xlim([min(days_vec) max(days_vec)]); ylim([0 3.0]);
ylabel('a.u.')
title(sprintf('P_{Z/G} in %cf=%.2f - %.2f Hz during %cT=%02d - %02d LT',char(916),f_rangeint(f_prec),f_rangeint(f_prec+1),char(916),LT_start,LT_end))

%% For thesis - direction (single date)

azim_plot=input('Enter a date in [yyyy,mm,dd] format ');
day_sel=find(floor(datenum(days_vec))==datenum(azim_plot));

f11=figure(11);

set(0,'DefaultFigureVisible','off');
histo1=polarhistogram(alpha3(:,day_sel),'BinEdges',binedge,'FaceColor','b');
azim_total=histo1.Values;

set(0,'DefaultFigureVisible','on');
if any(~isnan(alpha1(:,day_sel)))
    polarscatter(1,200,3,'o','MarkerFaceColor','m','MarkerEdgeColor','none');
    lgd=legend(sprintf('Upcoming EQs within %d days',daysafter_prec),'AutoUpdate','off','Location','southoutside');
    set(lgd,'FontSize',10);
    hold on
    for j=1:numel(azim_total)
        polarhistogram(mean(binedge(j:j+1))*ones(100,1),'BinEdges',binedge,...
            'FaceColor','k','FaceAlpha',(azim_total(j)-min(azim_total))/range(azim_total),'EdgeColor','none');
    end
    rlim([0 100]);
    for i=1:size(EQ_plrscat,1)
        if i==EQ_ID edge_alpha=1; else edge_alpha=0; end
        polarscatter(EQ_plrscat(i,1,day_sel),EQ_plrscat(i,4,day_sel),EQ_plrscat(i,3,day_sel),'o',...
            'MarkerFaceColor','m','MarkerEdgeColor','m',...
            'MarkerFaceAlpha',EQ_plrscat(i,5,day_sel),'MarkerEdgeAlpha',edge_alpha,'LineWidth',3);
    end
    hold off
    set(gca,'ThetaDir','clockwise','ThetaZeroLocation','top');
end
date_sel=char(days_vec(day_sel));
date_sel=date_sel(1:end-3);
title(sprintf('Azimuthal beam at %s on %s',stn,date_sel));

theta_actual=EQ_sel(EQ_ID,8);
disp(sprintf('Actual angle = %f',theta_actual));

%% 
load EQLIST_MANUSCRIPT
EQLIST_DATAOK=EQLIST_MANUSCRIPT;


for i=1:size(EQLIST_MANUSCRIPT,1)
    
  DATE_ALTERED0=EQLIST_DATAOK{i,1};
  DATE_ALTERED=datetime([str2num(strcat('20',DATE_ALTERED0(7:8))),str2num(DATE_ALTERED0(4:5)),str2num(DATE_ALTERED0(1:2))],'Format','dd/MM/yyyy');
  EQLIST_DATAOK{i,1}=char(DATE_ALTERED);
  
  LAT_ALTERED0=EQLIST_DATAOK{i,2};
  if LAT_ALTERED0<0 NORTHSOUTH='S'; else NORTHSOUTH='N'; end
  LAT_ALTERED1=abs(LAT_ALTERED0);
  LAT_ALTERED=sprintf('%06.2f%c %c',LAT_ALTERED1,char(176),NORTHSOUTH);
  EQLIST_DATAOK{i,2}=LAT_ALTERED;
    
  LON_ALTERED0=EQLIST_DATAOK{i,3};
  if LON_ALTERED0<0 EASTWEST='W'; else EASTWEST='E'; end
  LON_ALTERED1=abs(LON_ALTERED0);
  LON_ALTERED=sprintf('%06.2f%c %c',LON_ALTERED1,char(176),EASTWEST);
  EQLIST_DATAOK{i,3}=LON_ALTERED;
    
  MAG_ALTERED0=EQLIST_DATAOK{i,4};
  MAG_ALTERED=sprintf('%3.1f',MAG_ALTERED0);
  EQLIST_DATAOK{i,4}=MAG_ALTERED;
  
  DEP_ALTERED0=EQLIST_DATAOK{i,5};
  DEP_ALTERED=sprintf('%d',DEP_ALTERED0);
  EQLIST_DATAOK{i,5}=DEP_ALTERED;
  
  STN_ALTERED0=EQLIST_DATAOK{i,6};
  STN_ALTERED=sprintf('%s',STN_ALTERED0(1:3));
  EQLIST_DATAOK{i,6}=STN_ALTERED;
  
  DIS_ALTERED0=EQLIST_DATAOK{i,7};
  DIS_ALTERED=sprintf('%d',DIS_ALTERED0);
  EQLIST_DATAOK{i,7}=DIS_ALTERED;
  
  KLS_ALTERED0=EQLIST_DATAOK{i,8};
  KLS_ALTERED=sprintf('%d',round(KLS_ALTERED0,0));
  EQLIST_DATAOK{i,8}=KLS_ALTERED;
    
end

EQLIST_DATAOK2=EQLIST_DATAOK;

j=1;
for i=1:size(EQLIST_DATAOK2,1)
   if EQLIST_DATAOK2{j,9}==0
       EQLIST_DATAOK2(j,:)=[];
   else
       j=j+1;
   end
    
end



