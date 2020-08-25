%% Purpose
% Using polarization ellipse methodas in Ohta et. al (2013)
% Search for ULF precursor in frequency bands. Frequency might have been
% changed to higher ones (>0.1 Hz)
% Segmenting local time into 4 periods
% Include phase shift for determinancy (deltaS) calculation
%% Load initial variable
load VARIABLES_WORLD
open EQ_intrst
%% Place and date
EQ_num=50;                      %Select earthquake of interest #

stn0=EQ_intrst{EQ_num,1}; stn=stn0(1:3);
EQ_date0=EQ_intrst{EQ_num,3}; EQ_date=datenum(str2num(strcat('20',EQ_date0(7:8))),str2num(EQ_date0(4:5)),str2num(EQ_date0(1:2)));

date_start0=EQ_date-60;
date_end0=EQ_date+30;

date_start=datevec(date_start0);           
date_end=  datevec(date_end0);        
%% Customization
mag_min=5.0;                     %Minimum magnitude of earthquakes to be considered
dis_max=300;                     %Maximum epicentral distance from the station
depth_max=200;

LT_start=[22,00,02,04];
LT_end=[00,02,04,06];

f_rangeint=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10];
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
UT1m_sel=datevec(UT1m(86400*day_start-86399:86400*day_end));

H_dn=H(86400*day_start-86399:86400*day_end);
D_dn=D(86400*day_start-86399:86400*day_end);
Z_dn=Z(86400*day_start-86399:86400*day_end);
G_dn=sqrt(H_dn.^2+D_dn.^2);

%Removing noise/outlier
for i=1:5
    H_dn=medfilt1(H_dn,30,'omitnan','truncate');
    D_dn=medfilt1(D_dn,30,'omitnan','truncate');
    Z_dn=medfilt1(Z_dn,30,'omitnan','truncate');
    G_dn=medfilt1(G_dn,30,'omitnan','truncate');
end

for i=1:2
    
    H_sig=std(H_dn,'omitnan');
    H_mu=mean(H_dn,'omitnan');
    D_sig=std(D_dn,'omitnan');
    D_mu=mean(D_dn,'omitnan');
    Z_sig=std(Z_dn,'omitnan');
    Z_mu=mean(Z_dn,'omitnan');
    G_sig=std(G_dn,'omitnan');
    G_mu=mean(G_dn,'omitnan');
    
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
        if G_dn(j)>G_mu+5*G_sig||G_dn(j)<G_mu-5*G_sig
            G_dn(j)=NaN;
        end
    end
    
end

%% Nighttime data in LT

for l=1:numel(LT_start)
    LT_length(l)=LT_end(l)-LT_start(l);
    if LT_length(l)<0 LT_length(l)=LT_length(l)+24; end
    t_zone=timezone(stn_latlon(2),'degrees');
    t_start(l)=LT_start(l)+t_zone;
    if t_start(l)>24 t_start(l)=t_start(l)-24; end
    if t_start(l)<0 t_start(l)=t_start(l)+24;end
    t_starts(l)=t_start(l)*3600;
    t_int(l)=LT_length(l)*3600;
    
    for i=1:days_num
        H_night(:,i,l)=H_dn(t_starts(l)+1:t_starts(l)+t_int(l));
        D_night(:,i,l)=D_dn(t_starts(l)+1:t_starts(l)+t_int(l));
        Z_night(:,i,l)=Z_dn(t_starts(l)+1:t_starts(l)+t_int(l));
        G_night(:,i,l)=G_dn(t_starts(l)+1:t_starts(l)+t_int(l));
        UT1m_night(:,i,l)=UT1m(t_starts(l)+1:t_starts(l)+t_int(l));
        t_starts(l)=t_starts(l)+86400;
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

alpha1=NaN(size(H_night,1),days_num,numel(LT_start));
alpha2=NaN(size(H_night,1),days_num,numel(LT_start));

for l=1:numel(LT_start)
    for i=1:days_num
        
        % Precursor detection
        
        n_win=size(H_night,1);
        n_ovrlp=0;
        n_fft=n_win;
        fs=1;
        
        H_pass=filter1('bp',H_night(:,i,l),'fc',[f_rangeint(1) f_rangeint(end)],'fs',fs);
        D_pass=filter1('bp',D_night(:,i,l),'fc',[f_rangeint(1) f_rangeint(end)],'fs',fs);
        
        [Phh,fff]=periodogram(H_pass,hamming(n_win),n_fft,fs);
        [Pdd,fff]=periodogram(D_pass,hamming(n_win),n_fft,fs);
        
        [Phd,fff]=cpsd(H_pass,D_pass,hamming(n_win),n_ovrlp,n_fft,fs);
        [Pdh,fff]=cpsd(D_pass,H_pass,hamming(n_win),n_ovrlp,n_fft,fs);
        
        bin_floor=find(round(fff,4)==f_rangeint(1)); bin_ceil=find(round(fff,4)==f_rangeint(end));
        
        Phh1=Phh(bin_floor:bin_ceil,:);
        Pdd1=Pdd(bin_floor:bin_ceil,:);
        
        Phd1=Phd(bin_floor:bin_ceil,:);
        Pdh1=Pdh(bin_floor:bin_ceil,:);
        
        fff1=fff(bin_floor:bin_ceil,:);
        
      
        for a=1:numel(f_rangeint)
            f_rangeidx(a)=find(fff1==f_rangeint(a));
        end
        
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
        
        beta=0.5*asin(imag(Pdh_mean-Phd_mean) ./ ((Phh_mean-Pdd_mean).^2 + 4*Phh_mean.*Pdd_mean ).^0.5);
        
        deltaS(:,i,l)=real(((Phh_mean ./ Pdd_mean) - 1) ./ rms(tan(beta)));
        
        del_psi=30; psi=deg2rad(0:del_psi:180-del_psi);
        for f=1:numel(Phh_mean)
            Ptt=Phh_mean(f).*(cos(psi)).^2 + Pdd_mean(f).*(sin(psi)).^2 + real(Phd_mean(f)).*sin(2*psi);
            Prr=Phh_mean(f).*(sin(psi)).^2 + Pdd_mean(f).*(cos(psi)).^2 - real(Phd_mean(f)).*sin(2*psi);
     
            deltaSF_top=max(Ptt./Prr);
            deltaSF(f,i,l)=real((deltaSF_top - 1) / rms(tan(beta(f))));
        end
         
        % Direction estimation
        
        H_pass1=filter1('bp',H_night(:,i,l),'fc',[f_rangeint(1) f_rangeint(end)],'fs',fs);
        D_pass1=filter1('bp',D_night(:,i,l),'fc',[f_rangeint(1) f_rangeint(end)],'fs',fs);
        
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
        
        sum_HD=sqrt(H_pass.^2+D_pass.^2);
        mean_HD=sqrt(mean(H_pass.^2+D_pass.^2));
        
        for c=1:length(PE_alpha1)
            if sum_HD(c)>5*mean_HD
                alpha1(c,i,l)=PE_alpha1(c);
                alpha2(c,i,l)=PE_alpha2(c);
            end
        end
        
        
    end
end

alpha3=vertcat(alpha1,alpha2);
%% deltaS plotting

f1=figure(1);

stitle=sprintf('%s %s',stn,year);
suptitle(stitle)

sp1=subplot(5,1,1);
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

sp2=subplot(5,1,2);
surf(days_vec,f_rangeint(1:end-1),deltaS(:,:,1),'LineStyle','none')
view([0,0,90])
xlim([min(days_vec) max(days_vec)]);
ylim([f_rangeint(1) f_rangeint(end-1)]);
ylabel('f (Hz)')
title('\DeltaT = 22 - 00 LT')

sp3=subplot(5,1,3);
surf(days_vec,f_rangeint(1:end-1),deltaS(:,:,2),'LineStyle','none')
view([0,0,90])
xlim([min(days_vec) max(days_vec)]);
ylim([f_rangeint(1) f_rangeint(end-1)]);
ylabel('f (Hz)')
title('\DeltaT = 00 - 02 LT')

sp4=subplot(5,1,4);
surf(days_vec,f_rangeint(1:end-1),deltaS(:,:,3),'LineStyle','none')
view([0,0,90])
xlim([min(days_vec) max(days_vec)]);
ylim([f_rangeint(1) f_rangeint(end-1)]);
ylabel('f (Hz)')
title('\DeltaT = 02 - 04 LT')

sp5=subplot(5,1,5);
surf(days_vec,f_rangeint(1:end-1),deltaS(:,:,4),'LineStyle','none')
view([0,0,90])
xlim([min(days_vec) max(days_vec)]);
ylim([f_rangeint(1) f_rangeint(end-1)]);
ylabel('f (Hz)')
title('\DeltaT = 04 - 06 LT')

% colmap=gray;
% colmap=flipud(colmap);
% colormap(colmap);

linkaxes([sp1,sp2,sp3,sp4,sp5],'x')
%% Azimuthal angle polar plotting

azim_plot=input('Enter a date in [yyyy,mm,dd] format ');

f2=figure(2);

day_sel=find(datenum(days_vec)==datenum(azim_plot));

subplot(2,2,1)
if any(~isnan(alpha3(:,day_sel,1)))
    hist1=polarhistogram(alpha3(:,day_sel,1),36,'FaceColor','b');
    hold on
    polarscatter(deg2rad(EQ_sel(:,8)),max(hist1.Values)*ones(1,size(EQ_sel,1)),'k^','MarkerFaceColor','y')
    polarscatter(deg2rad(max(EQ_sel(:,8))),max(hist1.Values),'k^','MarkerFaceColor','r')
    hold off
    pax = gca;
    pax.ThetaDir = 'clockwise';
    pax.ThetaZeroLocation = 'top';
end

subplot(2,2,2)
if any(~isnan(alpha3(:,day_sel,2)))
    hist2=polarhistogram(alpha3(:,day_sel,2),36,'FaceColor','b');
    hold on
    polarscatter(deg2rad(EQ_sel(:,8)),max(hist2.Values)*ones(1,size(EQ_sel,1)),'k^','MarkerFaceColor','y')
    hold off
    pax = gca;
    pax.ThetaDir = 'clockwise';
    pax.ThetaZeroLocation = 'top';
end

subplot(2,2,3)
if any(~isnan(alpha3(:,day_sel,3)))
    hist3=polarhistogram(alpha3(:,day_sel,3),36,'FaceColor','b');
    hold on
    polarscatter(deg2rad(EQ_sel(:,8)),max(hist3.Values)*ones(1,size(EQ_sel,1)),'k^','MarkerFaceColor','y')
    hold off
    pax = gca;
    pax.ThetaDir = 'clockwise';
    pax.ThetaZeroLocation = 'top';
end

subplot(2,2,4)
if any(~isnan(alpha3(:,day_sel,4)))
    hist4=polarhistogram(alpha3(:,day_sel,4),36,'FaceColor','b');
    hold on
    polarscatter(deg2rad(EQ_sel(:,8)),max(hist4.Values)*ones(1,size(EQ_sel,1)),'k^','MarkerFaceColor','y')
    hold off
    pax = gca;
    pax.ThetaDir = 'clockwise';
    pax.ThetaZeroLocation = 'top';
end



%% Raw data plotting

f3=figure(3)
sec_vec=datetime(date_start1:1/86400:date_end1,'ConvertFrom','datenum','Format','dd/MM/yyyy');
sp1=subplot(3,1,1)
plot(sec_vec,H_dn)
ylabel('H (nT)')
sp2=subplot(3,1,2)
plot(sec_vec,D_dn)
ylabel('D (nT)')
sp3=subplot(3,1,3)
plot(hourly_vec,ap)
hold on
plot(hourly_vec,Dst)
ylabel('ap   Dst')
xlim([min(hourly_vec) max(hourly_vec)])
hold off
linkaxes([sp1,sp2,sp3],'x')







