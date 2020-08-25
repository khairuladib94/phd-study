%% %Purpose
%To plot daily Z/G and geomagnetic indices temporal evolution with
%customizable segmentation period (hourly as default) and nighttime period labels
%No exponential fitting
%Perform Welch in all four frequency channels at a same time
%Calculate azimuthal angle and amplitude of anomalies by assigning
%different colors for different earthquakes
%OPTIONAL:
%-Inspect raw data
%% %Place and time inputs
stn='TNO';                         %Abbreviation of station name
date_start=[2010,11,01];           %Insert custom start and end dates
date_end=  [2011,05,01];           %Period spanning through 3 consecutive years is the maximum
%% %Customization
mag_min=5.0;                     %Minimum magnitude of earthquakes to be considered
dis_max=180;                     %Maximum epicentral distance from the station
depth_max=200;                   %Maximum hypocentral depth

%Frequency channels
f_1=0.0067;
f_2=0.0100;

f_3=0.0100;
f_4=0.0220;

f_5=0.0220;
f_6=0.0500;

f_7=0.0500;
f_8=0.1000;

N_seg=3600;                     %Segmentation length
pfit_order=0;                   %Order of polynomial fit for ZG detrending
movmean_val=86400/N_seg;        %Moving mean of ZG and indices

%Normalization type
norm_type='mo';                    %'mo'=monthly, 'eq'=equal or 'no'=no normalization

%Nighttime indicator patch
LT_start=22;
LT_end=04;
%% %Time period

date_start1=datenum(horzcat(date_start,[0,0,0]));
date_end1=datenum(horzcat(date_end,[23,59,59]));
days_vec=datetime(datevec(date_start1:1/(86400/N_seg):date_end1),'Format','dd/MM/yyyy HH');
date_vec=datevec(date_start1:1/(86400/N_seg):date_end1);
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
        if EQ_angle>180
            EQ_angle=EQ_angle-360;
        end
        
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
%% %Geomagnetic data pre-processing
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
for i=1:3
    H_dn=medfilt1(H_dn,50,'omitnan','truncate');
    D_dn=medfilt1(D_dn,50,'omitnan','truncate');
    Z_dn=medfilt1(Z_dn,50,'omitnan','truncate');
    G_dn=medfilt1(G_dn,50,'omitnan','truncate');
end

for i=1:5
    
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

if norm_type=='mo'
    %Normalization based on months
    for i=1:size(year_vec,2)
        same_yr=find(UT1m_sel(:,1)==year_vec(i));
        mon_vec=unique(UT1m_sel(same_yr,2),'stable');
        for j=1:size(mon_vec,1)
            same_mon=find(UT1m_sel(:,2)==mon_vec(j));
            H_dn1(same_mon,1)=normalize(H_dn(same_mon,1));
            D_dn1(same_mon,1)=normalize(D_dn(same_mon,1));
            G_dn1(same_mon,1)=normalize(G_dn(same_mon,1));
            Z_dn1(same_mon,1)=normalize(Z_dn(same_mon,1));
        end
    end
elseif norm_type=='eq'
    %Equal normalization
    H_dn1=normalize(H_dn);
    D_dn1=normalize(D_dn);
    Z_dn1=normalize(Z_dn);
    G_dn1=normalize(G_dn);
elseif norm_type=='no'
    %No normalization
    H_dn1=H_dn;
    D_dn1=D_dn;
    Z_dn1=Z_dn;
    G_dn1=G_dn;
end

H_seg=reshape(H_dn1,N_seg,[]);
D_seg=reshape(D_dn1,N_seg,[]);
Z_seg=reshape(Z_dn1,N_seg,[]);
G_seg=reshape(G_dn1,N_seg,[]);

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
        AP(k,1)=index_geomag(i,6);
        DST(k,1)=index_geomag(i,5);
        k=k+1;
    end
end

idx_num=N_seg/3600;
if idx_num>1
    ap=reshape(AP,idx_num,[]);
    Dst=reshape(DST,idx_num,[]);
    ap=mean(ap)';
    Dst=mean(Dst)';
elseif idx_num==1
    ap=AP;
    Dst=DST;
elseif idx_num<1
    idx_mult=days_num/length(AP);
    j=1;
    for i=1:length(AP)
        ap(j:j+idx_mult-1,1)=AP(i,1);
        Dst(j:j+idx_mult-1,1)=DST(i,1);
        j=j+idx_mult;
    end
end
%% %PSD calculation
n_win=1800;
n_ovrlp=0.5*n_win;
n_fft=n_win;
fs=1;

f_bot=[f_1,f_3,f_5,f_7];
f_top=[f_2,f_4,f_6,f_8];

for i=1:4
    bin_1=[];
    bin_2=[];
    
    [H_wlch,f_wlch]=pwelch(H_seg,hamming(n_win),n_ovrlp,n_fft,fs);
    [D_wlch,f_wlch]=pwelch(D_seg,hamming(n_win),n_ovrlp,n_fft,fs);
    [Z_wlch,f_wlch]=pwelch(Z_seg,hamming(n_win),n_ovrlp,n_fft,fs);
    [G_wlch,f_wlch]=pwelch(G_seg,hamming(n_win),n_ovrlp,n_fft,fs);
    
    f_wlch4=round(f_wlch,4);
    bin_1=find(f_wlch4==f_bot(i));
    bin_2=find(f_wlch4==f_top(i));
    if isempty(bin_1) || isempty(bin_2)
        f_wlch3=round(f_wlch,3);
        bin_1=floor(mean(find(f_wlch3==f_bot(i))));
        bin_2=floor(mean(find(f_wlch3==f_top(i))));
    end
    
    H_wlchf=H_wlch(bin_1:bin_2,:);
    H_wlchmu(:,:,i)=mean(H_wlchf,1,'omitnan');
    D_wlchf=D_wlch(bin_1:bin_2,:);
    D_wlchmu(:,:,i)=mean(D_wlchf,1,'omitnan');
    Z_wlchf=Z_wlch(bin_1:bin_2,:);
    Z_wlchmu(:,:,i)=mean(Z_wlchf,1,'omitnan');
    G_wlchf=G_wlch(bin_1:bin_2,:);
    G_wlchmu(:,:,i)=mean(G_wlchf,1,'omitnan');
end

H_norm=H_wlchmu;
D_norm=D_wlchmu;
Z_norm=Z_wlchmu;
G_norm=G_wlchmu;

G_norm(G_norm==0)=NaN;
%% %Power ratio

ZG=Z_norm./G_norm;

for i=1:4
    ZG_fit=polyfit(1:days_num,fillgaps(ZG(:,:,i)),pfit_order);
    ZG_curve(:,:,i)=polyval(ZG_fit,1:days_num);
    if pfit_order>0 ZG(:,:,i)=ZG(:,:,i)-ZG_curve(:,:,i); end
end

%mu+/-2sigma
ZG_mu=nanmean(ZG);
ZG_sig=nanstd(ZG);
ZG_mup2sig=ones(1,days_num).*(ZG_mu+2*ZG_sig);
ZG_mum2sig=ones(1,days_num).*(ZG_mu-2*ZG_sig);
%% %Disturbed days and anomalies detection
dstb=NaN(1,days_num,4);
anom=NaN(1,days_num,4);
for j=1:4
    for i=1:days_num
        if ZG(1,i,j)>ZG_mup2sig(1,1,j) || ZG(1,i,j)<ZG_mum2sig(1,1,j)
            if ap(i,1)>50 || Dst(i,1)<-50
                dstb(1,i,j)=ZG(1,i,j);
                if i<days_num
                    dstb(1,i+1,j)=ZG(1,i+1,j);
                end
                if i>1
                    dstb(1,i-1,j)=ZG(1,i-1,j);
                end
            else
                anom(1,i,j)=ZG(1,i,j);
                if i<days_num
                    anom(1,i+1,j)=ZG(1,i+1,j);
                end
                if i>1
                    anom(1,i-1,j)=ZG(1,i-1,j);
                end
            end
        end
    end
end

largeEQ=EQ_sel(find(EQ_sel(:,4)==max(EQ_sel(:,4))),1);
largeEQ_3D=datetime(largeEQ-3.0,'ConvertFrom','datenum','Format','dd/MM/yyyy HH');
largeEQ_1W=datetime(largeEQ-7.0,'ConvertFrom','datenum','Format','dd/MM/yyyy HH');
largeEQ_2W=datetime(largeEQ-14.0,'ConvertFrom','datenum','Format','dd/MM/yyyy HH');
largeEQ_1M=datetime(largeEQ-30.0,'ConvertFrom','datenum','Format','dd/MM/yyyy HH');
%% %Azimuthal angle calculation
% XX=H_seg;
% YY=D_seg;
% ZZ=Z_seg;
% n_win_dir=120;
% n_ovrlp_dir=60;
% n_fft_dir=n_win_dir;
% for i=1:size(XX,2)
%     XX_spec=spectrogram(XX(:,i),'yaxis',hamming(n_win_dir),n_ovrlp_dir,n_fft_dir,1);
%     YY_spec=spectrogram(YY(:,i),'yaxis',hamming(n_win_dir),n_ovrlp_dir,n_fft_dir,1);
%     ZZ_spec=spectrogram(ZZ(:,i),'yaxis',hamming(n_win_dir),n_ovrlp_dir,n_fft_dir,1);
%     bin036=round(0.36/(0.5/size(XX_spec,1)));
%     XX_036=real(XX_spec(bin036,:));
%     YY_036=real(YY_spec(bin036,:));
%     ZZ_036=real(ZZ_spec(bin036,:));
%     XY_mat=[XX_036;YY_036];
%     ZZ_mat=(ZZ_036');
%     AB=real((inv(XY_mat*XY_mat'))*(XY_mat*ZZ_mat));
%     AB_amp(1,i)=sqrt(AB(1)^2+AB(2)^2);
%     AB_theta(1,i)=atan2d(AB(2),AB(1));
% end
XX=H_wlchf;
YY=D_wlchf;
ZZ=Z_wlchf;
for i=1:size(XX,2)
    XY_mat=[XX(:,i) YY(:,i)]';
    ZZ_mat=ZZ(:,i);
    AB=real((inv(XY_mat*XY_mat'))*(XY_mat*ZZ_mat));
    AB_amp(1,i)=(sqrt(AB(1)^2+AB(2)^2));
    AB_theta(1,i)=atan2d(AB(2),AB(1));
end
%% %Figure 1
title_sup=sprintf('%s station, %s - %s',station,datetime(datevec(datenum(date_start)),'Format','dd/MM/yyyy'),datetime(datevec(datenum(date_end)),'Format','dd/MM/yyyy'));

f1=figure(2);
stitle=suptitle(sprintf('%s',title_sup));
stitle_pos =get(stitle,'position');
stitle_pos=[stitle_pos(1) stitle_pos(2)+0.01 stitle_pos(3)];
set(stitle,'position',stitle_pos);

sp1=subplot(3,1,1)
plot(days_vec,ap,'k');
hold on
plot(days_vec,Dst,'r');
plot(days_vec,ones(length(days_vec),1)*(50),'k--');
plot(days_vec,ones(length(days_vec),1)*(-50),'r--');
hold off
ylabel('Dst               ap');
if max(ap)>50 && min(Dst)<-50 ylim([min(Dst)-10 max(ap)+10])
elseif max(ap)>50 ylim([-60 max(ap)+10])
elseif min(Dst)<-50 ylim([min(Dst)-10 60])
else ylim([-60 60])
end
yyaxis right
set(gca, 'YColor', 'b')
for i=1:size(EQ_sel,1)
    hold on
    j=datetime(datevec(EQ_sel(i,1)),'Format','dd/MM/yyyy HH:mm:ss');
    if 0.025*EQ_sel(i,2)<EQ_sel(i,3)-4.5
        plot(j,EQ_sel(i,4),'^','MarkerEdgeColor',EQ_sel(i,9:11),'MarkerFaceColor',EQ_sel(i,9:11),'MarkerSize',5);
        label={sprintf('%s',string(datetime(j,'Format','dd/MM/yyyy')));
            sprintf('M%0.1f D=\t%.0fKM',EQ_sel(i,3),EQ_sel(i,5));
            sprintf('d=%.0fKM',EQ_sel(i,2))};
        text(j,EQ_sel(i,4),label,'VerticalAlignment','top','HorizontalAlignment','left','FontSize',9)
    else
        plot(j,EQ_sel(i,4),'o','MarkerEdgeColor',EQ_sel(i,9:11),'MarkerFaceColor',EQ_sel(i,9:11),'MarkerSize',5);
    end
    if EQ_sel(i,4)==max(EQ_sel(:,4))
        u=j;
        t=i;
        plot([u u],[-1e5 1e5],'-','Color',EQ_sel(t,9:11));
        label={sprintf('%s',string(datetime(j,'Format','dd/MM/yyyy')));
            sprintf('M%0.1f D=\t%.0fKM',EQ_sel(i,3),EQ_sel(i,5));
            sprintf('d=%.0fKM',EQ_sel(i,2))};
        text(j,EQ_sel(i,4),label,'VerticalAlignment','top','HorizontalAlignment','left','FontSize',9)
    end
end
hold off
xlim([min(days_vec) max(days_vec)])
if ~any(isnan(EQ_sel)) ylim([0 1.1*max(EQ_sel(:,4))]); end
ylabel('K_{LS}');
title('Geomagnetic (left) & local seismicity (right) indices');

sp2=subplot(3,1,2)
hold on
for i=1:4
    plot(days_vec,ZG(:,:,i));
end
legend({sprintf('\\Deltaf=%.3f-%.3f Hz',f_bot(1),f_top(1)),sprintf('\\Deltaf=%.3f-%.3f Hz',f_bot(2),f_top(2)),sprintf('\\Deltaf=%.3f-%.3f Hz',f_bot(3),f_top(3)),sprintf('\\Deltaf=%.3f-%.3f Hz',f_bot(4),f_top(4))},'location','northwestoutside','orientation','horizontal');
for i=1:3
    plot(days_vec,anom(:,:,i),'r','HandleVisibility','off');
    plot(days_vec,dstb(:,:,i),'g','HandleVisibility','off');
end

plot(days_vec,anom(:,:,4),'r');
plot(days_vec,dstb(:,:,4),'g');

if ~any(isnan(EQ_sel))
    plot([u u],[-1e5 1e5],'-','Color',EQ_sel(t,9:11),'HandleVisibility','off');
    plot([largeEQ_3D largeEQ_3D],[-1e5 1e5],'--','Color',[0.5 0.5 0.5],'HandleVisibility','off');
    text(largeEQ_3D,min(min(ZG)),'-3D','VerticalAlignment','top','HorizontalAlignment','right','FontSize',9)
    plot([largeEQ_1W largeEQ_1W],[-1e5 1e5],'--','Color',[0.5 0.5 0.5],'HandleVisibility','off');
    text(largeEQ_1W,min(min(ZG)),'-1W','VerticalAlignment','top','HorizontalAlignment','right','FontSize',9)
    plot([largeEQ_2W largeEQ_2W],[-1e5 1e5],'--','Color',[0.5 0.5 0.5],'HandleVisibility','off');
    text(largeEQ_2W,min(min(ZG)),'-2W','VerticalAlignment','top','HorizontalAlignment','right','FontSize',9)
    plot([largeEQ_1M largeEQ_1M],[-1e5 1e5],'--','Color',[0.5 0.5 0.5],'HandleVisibility','off');
    text(largeEQ_1M,min(min(ZG)),'-1M','VerticalAlignment','top','HorizontalAlignment','right','FontSize',9)
end
hold off
legend({sprintf('\\Deltaf=%.3f-%.3f Hz',f_bot(1),f_top(1)),sprintf('\\Deltaf=%.3f-%.3f Hz',f_bot(2),f_top(2)),sprintf('\\Deltaf=%.3f-%.3f Hz',f_bot(3),f_top(3)),sprintf('\\Deltaf=%.3f-%.3f Hz',f_bot(4),f_top(4)),'Anomaly','Disturbed',},'location','northeast','orientation','horizontal');
ylabel('Z/G power ratio');
if 0.8*min(min(ZG))<0 ylim([0.8*min(min(ZG)) 1.2*max(max(ZG))]);
else ylim([0 1.2*max(max(ZG))]); end
xlim([min(days_vec) max(days_vec)])
title('Z/G polarisation ratio');

sp3=subplot(3,1,3)
hold on
plotcol={'r','b','g','c'};
for i=1:4
    plot(days_vec,Z_norm(:,:,i),plotcol{i},'LineStyle','-');
end
ylim([0 0.5]); ylabel('Z_{power} (upright)')
yyaxis right
set(gca,'YDir','reverse');
for i=1:4
    plot(days_vec,G_norm(:,:,i),plotcol{i},'LineStyle','-','HandleVisibility','off');
end
ylim([0 0.5]); ylabel('G_{power} (inverted)')
legend({sprintf('\\Deltaf=%.3f-%.3f Hz',f_bot(1),f_top(1)),sprintf('\\Deltaf=%.3f-%.3f Hz',f_bot(2),f_top(2)),sprintf('\\Deltaf=%.3f-%.3f Hz',f_bot(3),f_top(3)),sprintf('\\Deltaf=%.3f-%.3f Hz',f_bot(4),f_top(4))},'location','northeast','orientation','horizontal');
hold off


% sp5=subplot(5,1,5)
% plot(days_vec,movmean(AB_theta,1),'k')
% hold on
% for i=1:size(EQ_sel,1)
%     j=datetime(datevec(EQ_sel(i,1)),'Format','dd/MM/yyyy HH:mm:ss');
%     plot(j,EQ_sel(i,8),'o','MarkerEdgeColor',EQ_sel(i,9:11),'MarkerFaceColor',EQ_sel(i,9:11),'MarkerSize',5);
%     if EQ_sel(i,4)==max(EQ_sel(:,4))
%         plot(days_vec,ones(1,days_num)*EQ_sel(i,8),'--','Color',EQ_sel(i,9:11))
%     end
% end
% ylabel('Azimuthal angle (\theta)');
% hold off
%
% yyaxis right
% plot(days_vec,AB_amp)
% ylabel('Amplitude');
% xlim([min(days_vec) max(days_vec)])
% title('Source azimuthal angle (left) and amplitude (right)');

descr=sprintf('Maximum distance = %d km | Minimum magnitude = M%.1f | Maximum depth = %d km | Segmentation = %d s | Moving mean = %d | Detrending polynomial order = %d',dis_max,mag_min,depth_max,N_seg,movmean_val,pfit_order);
text('String',descr,'Position',[0.5 1.12],'Units','normalized','EdgeColor','k','HorizontalAlignment','center');

linkaxes([sp1,sp2,sp3],'x')

set(gcf, 'Position', get(0, 'Screensize'));
%% %Figure 2
% f2=figure(2);
%
% worldmap([stn_latlon(1)-8 stn_latlon(1)+8],[stn_latlon(2)-8 stn_latlon(2)+8])
% geoshow('landareas.shp','FaceColor','White')
% title(sprintf('Earthquake Map | %s - %s, M > %0.1f, d < %.0f km',datetime(datevec(datenum(date_start)),'Format','dd/MM/yyyy'),datetime(datevec(datenum(date_end)),'Format','dd/MM/yyyy'),mag_min,dis_max));
%
% if ~isnan(mean(EQ_sel))
%     for i=1:size(EQ_sel,1)
%         EQ_lat=EQ_sel(i,6);
%         EQ_lon=EQ_sel(i,7);
%         EQ_radius=((EQ_sel(i,3)-4.5)/0.025);
%         date_text=datetime(datevec(EQ_sel(i,1)),'Format','dd/MM/yyyy');
%         circlem(EQ_lat,EQ_lon,EQ_radius,'edgecolor','none','facecolor',colour,'facealpha','0.2');
%         if (EQ_sel(i,3)>=6.0 && EQ_sel(i,2)<=100) || EQ_sel(i
%         if EQ_sel(i,5)<=50
%             colour='r';
%         elseif (EQ_sel(i,5)>50) && (EQ_sel(i,5)<=200)
%             colour='b';
%         else EQ_sel(i,5)>200
%             colour='g';
%         end,2)<=30 || EQ_sel(i,3)>=7.0
%             textm(EQ_lat,EQ_lon,sprintf('%s(%0.1f)',date_text,EQ_sel(i,3)),'FontSize',8);
%         end
%     end
% end
% plotm(stn_latlon(1),stn_latlon(2),'k^','MarkerFaceColor','k','MarkerSize',5);
% textm(stn_latlon(1),stn_latlon(2)+0.2,stn);
% for i=1:length(stn_MAGDAS)
%     if stn_MAGDAS(i,2)>stn_latlon(1)-8 && stn_MAGDAS(i,2)<stn_latlon(1)+8 && stn_MAGDAS(i,3)>stn_latlon(2)-8 && stn_MAGDAS(i,3)<stn_latlon(2)+8
%         stnnearby_latlon=stn_MAGDAS(i,2:3);
%         stnnearby_name=stn_vec(i);
%         plotm(stnnearby_latlon(1),stnnearby_latlon(2),'k^','MarkerSize',5);
%         textm(stnnearby_latlon(1),stnnearby_latlon(2)+0.2,stnnearby_name);
%     end
% end
%
% set(gcf, 'Position', get(0, 'Screensize'));
%% %Saving figures 1&2
ddMMyy1=upper(string(datetime(min(days_vec),'Format','ddMMyy')));
ddMMyy2=upper(string(datetime(max(days_vec),'Format','ddMMyy')));

today=char(datetime('today','Format','dd-MM-yyyy'));
today=strcat(today,'\Adib_10\');
path1='D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\';
mkdir(fullfile(path1,today));

path2=strcat('D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\',today);

figname1=char(strcat(stn,'_',ddMMyy1,'-',ddMMyy2,'(ZG)'));
figname2=char(strcat(stn,'_',ddMMyy1,'-',ddMMyy2,'(MAP)'));

figcheck=exist(strcat(path2,figname1,'.tif'));
if figcheck==2
    askoverwrite=input('Figure exists. Overwrite? Enter 0=Don''t save, 1=Overwrite or 2=Save with different name ');
    if askoverwrite==0
        disp('Figure not saved')
    elseif askoverwrite==1 saveas(f1,fullfile(path2,figname1),'tiff');
        disp('Figure saved')
    elseif askoverwrite==2
        for i=1:100
            if exist(strcat(path2,figname1,sprintf('-%d',i),'.tif'))~=2
                figname1ext=strcat(figname1,char(sprintf('-%d',i)));
                saveas(f1,fullfile(path2,figname1ext),'tiff');
                disp('Figure saved')
                break
            end
        end
    end
else
    asksave=input('Save figure? 0=Don''t save and 1=Save ');
    if asksave==0
        disp('Figure not saved')
    elseif asksave==1
        saveas(f1,fullfile(path2,figname1),'tiff'); disp('Figure saved')
    end
end
% saveas(f2,fullfile(path2,figname2),'tiff');
%% %Inspect data

%Input
inspect=string(sprintf('%08d',input('Inspect data? Enter 0=No, 1=Yes, similar period or DDMMYYYY=Yes, more specific. Start date? ')));
if inspect=="00000000"
    disp('Analysis ends')
elseif inspect=="00000001"
    datenum_instart=date_start1;
    datenum_insend=date_end1;
else
    date_instart=inspect;
    inspect2=string(sprintf('%08d',input('End date? ')));
    if inspect2=="00000001"
        date_insend=date_instart;
        datenum_instart=datenum(date_instart,'ddmmyyyy');
        datenum_insend=datenum(strcat(date_instart,'235959'),'ddmmyyyyhhMMss');
    else
        date_insend=inspect2;
        datenum_instart=datenum(date_instart,'ddmmyyyy');
        datenum_insend=datenum(strcat(date_insend,'235959'),'ddmmyyyyhhMMss');
    end
    
end

%Proceed analysis
if inspect~="00000000"
    
    %Data selection
    instart=datenum_instart-data_start+1;
    insend=floor(datenum_insend)-data_start+1;
    insdatenum_vec=datenum_instart:datenum_insend;
    insdays_vec=datetime(datevec(datenum_instart:1/86400:datenum_insend),'Format','dd/MM/yyyy HH:mm:ss');
    
    G=sqrt(H.^2+D.^2);
    G_ins=G(86400*instart-86399:86400*insend);
    Z_ins=Z(86400*instart-86399:86400*insend);
    
    %Removing noise/outlier
    for i=1:3
        G_ins=medfilt1(G_ins,50,'omitnan','truncate');
        Z_ins=medfilt1(Z_ins,50,'omitnan','truncate');
    end
    
    for i=1:5
        
        G_sig=std(G_ins,'omitnan');
        G_mu=mean(G_ins,'omitnan');
        Z_sig=std(Z_ins,'omitnan');
        Z_mu=mean(Z_ins,'omitnan');
        
        for j=1:length(G_ins)
            
            if G_ins(j)>G_mu+5*G_sig||G_ins(j)<G_mu-5*G_sig
                G_ins(j)=NaN;
            end
            if Z_ins(j)>Z_mu+5*Z_sig||Z_ins(j)<Z_mu-5*Z_sig
                Z_ins(j)=NaN;
            end
            
        end
        
    end
    
    %Geomagnetic indices
    clearvars AP_ins ap_ins DST_ins Dst_ins
    k=1;
    for i=1:length(index_geomag)
        m=datenum(index_geomag(i,1),1,index_geomag(i,2));
        if m>=datenum_instart && m<=datenum_insend
            AP_ins(k,1)=index_geomag(i,6);
            DST_ins(k,1)=index_geomag(i,5);
            k=k+1;
        end
    end
    
    j=1;
    for i=1:length(AP_ins)
        ap_ins(j:j+3600-1,1)=AP_ins(i,1);
        Dst_ins(j:j+3600-1,1)=DST_ins(i,1);
        j=j+3600;
    end
    
    %Diff calculation
    G_diff=diff(G_ins);
    G_diff(length(insdays_vec))=NaN;
    Z_diff=diff(Z_ins);
    Z_diff(length(insdays_vec))=NaN;
    
    %Spectrogram parameters
    n_win=1800;
    n_ovrlp=0.5*n_win;
    n_fft=n_win;
    fs=1;
    
    %Plotting
    title_sup=sprintf('Data inspection at %s station, %s - %s',station,datetime(datevec(datenum_instart),'Format','dd/MM/yyyy'),datetime(datevec(floor(datenum_insend)),'Format','dd/MM/yyyy'));
    
    f3=figure(3);
    stitle=suptitle(sprintf('%s',title_sup));
    stitle_pos =get(stitle,'position');
    stitle_pos=[stitle_pos(1) stitle_pos(2)+0.01 stitle_pos(3)];
    set(stitle,'position',stitle_pos);
    
    spA=subplot(3,2,1);
    plot(insdays_vec,G_ins);
    ylabel('G_{raw}')
    xlim([min(insdays_vec) max(insdays_vec)])
    yyaxis right
    set(gca, 'YColor', 'b')
    for i=1:size(EQ_sel,1)
        hold on
        j=datetime(datevec(EQ_sel(i,1)),'Format','dd/MM/yyyy HH:mm:ss');
        if 0.025*EQ_sel(i,2)<EQ_sel(i,3)-4.5
            plot(j,EQ_sel(i,4),'^','MarkerEdgeColor',EQ_sel(i,9:11),'MarkerFaceColor',EQ_sel(i,9:11),'MarkerSize',5);
        else
            plot(j,EQ_sel(i,4),'o','MarkerEdgeColor',EQ_sel(i,9:11),'MarkerFaceColor',EQ_sel(i,9:11),'MarkerSize',5);
            
        end
    end
    hold off
    title('Raw G component');
    
    spB=subplot(3,2,2);
    plot(insdays_vec,Z_ins);
    ylabel('Z_{raw}')
    xlim([min(insdays_vec) max(insdays_vec)])
    yyaxis right
    set(gca, 'YColor', 'b')
    for i=1:size(EQ_sel,1)
        hold on
        j=datetime(datevec(EQ_sel(i,1)),'Format','dd/MM/yyyy HH:mm:ss');
        if 0.025*EQ_sel(i,2)<EQ_sel(i,3)-4.5
            plot(j,EQ_sel(i,4),'^','MarkerEdgeColor',EQ_sel(i,9:11),'MarkerFaceColor',EQ_sel(i,9:11),'MarkerSize',5);
        else
            plot(j,EQ_sel(i,4),'o','MarkerEdgeColor',EQ_sel(i,9:11),'MarkerFaceColor',EQ_sel(i,9:11),'MarkerSize',5);
        end
    end
    hold off
    title('Raw Z component');
    
    spC=subplot(3,2,3);
    plot(insdays_vec,G_diff);
    ylabel('G_{diff}')
    ylim([-3 3])
    xlim([min(insdays_vec) max(insdays_vec)])
    yyaxis right
    plot(insdays_vec,Dst_ins);
    ylim([-inf inf])
    ylabel('Dst')
    title('Diff(G) & Dst index');
    
    spD=subplot(3,2,4);
    plot(insdays_vec,Z_diff);
    ylabel('Z_{diff}')
    ylim([-3 3])
    xlim([min(insdays_vec) max(insdays_vec)])
    yyaxis right
    plot(insdays_vec,ap_ins);
    ylim([-inf inf])
    ylabel('ap')
    title('Diff(Z) & ap index');
    
    spE=subplot(3,2,5);
    [Sspec,Fspec,Tspec]=spectrogram(G_ins,'yaxis',hamming(n_win),n_ovrlp,n_fft,1);
    Tspec_new=insdays_vec(Tspec,1);
    surf(Tspec_new,Fspec,abs(Sspec),'EdgeColor','none');
    axis xy; axis tight; colormap(jet); view(0,90);
    ylim([0 0.1]); xlim([Tspec_new(1) Tspec_new(end)])
    ylabel('Frequency (Hz)')
    caxis([0 30])
    colorbar('off')
    colorbar('east')
    title('Spectrogram of G');
    
    spF=subplot(3,2,6);
    [Sspec,Fspec,Tspec]=spectrogram(Z_ins,'yaxis',hamming(n_win),n_ovrlp,n_fft,1);
    Tspec_new=insdays_vec(Tspec,1);
    surf(Tspec_new,Fspec,abs(Sspec),'EdgeColor','none');
    axis xy; axis tight; colormap(jet); view(0,90);
    ylim([0 0.1]); xlim([Tspec_new(1) Tspec_new(end)])
    ylabel('Frequency (Hz)')
    caxis([0 30])
    colorbar('off')
    colorbar('east')
    title('Spectrogram of Z');
    
    linkaxes([spA,spB,spC,spD,spE,spF],'x')
    linkaxes([spC,spD],'y')
    linkaxes([spE,spF],'y')
    set(gcf, 'Position', get(0, 'Screensize'));
    
    %Saving figure
    insddMMyy1=upper(string(datetime(min(insdays_vec),'Format','ddMMyy')));
    insddMMyy2=upper(string(datetime(max(insdays_vec),'Format','ddMMyy')));
    
    today=char(datetime('today','Format','dd-MM-yyyy'));
    today=strcat(today,'\Adib_10\');
    path1='D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\';
    mkdir(fullfile(path1,today));
    
    path2=strcat('D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\',today);
    
    figname3=char(strcat(stn,'_',insddMMyy1,'-',insddMMyy2,'(INSPECT)'));
    
    asksave=input('Save figure? 0=Don''t save and 1=Save ');
    if asksave==0
        disp('Figure not saved')
    elseif asksave==1
        saveas(f3,fullfile(path2,figname3),'tiff');
        disp('Figure saved')
    end
    
    disp('Analysis ends');
end
