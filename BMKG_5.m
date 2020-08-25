%% %Purpose
%To inspect raw data, plot maximum Z/H plot in each window and finally plot
%bandpassed signal
%% %Place and time inputs
stn='ONW';                                 %Abbreviation of station name
date_start=[2008,05,15,00,00,00];           %Insert custom start and end dates
date_end=  [2008,06,15,23,59,59];           %Period spanning through 3 consecutive years is the maximum  
addtitle='';
%% %Customization
%Set magnitude % distance
mag_min=5.0;                     
dis_max=300;    

%Segmentation (in secs)
N_seg=3600;

%Frequency range for PSD
f_1=0.05;
f_2=0.1;

%Night time setup
LT_start=0;
LT_end=6;
%% %Time vector
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

%Creating time vector
time_vec=datetime(datevec(datenum_start:1/86400:datenum_end),'Format','dd/MM/yyyy HH:mm:ss');
timenum_vec=datenum_start:1/86400:datenum_end;
welchnum_vec=datetime(datevec(datenum_start:1/(86400/N_seg):datenum_end),'Format','dd/MM/yyyy HH');
datenum_vec=datenum_start:1:datenum_end;
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
G_period=sqrt(H_period.^2+D_period.^2);

%Removing noise/outlier
H_dn=medfilt1(H_period,'omitnan');
D_dn=medfilt1(D_period,'omitnan');
Z_dn=medfilt1(Z_period,'omitnan');
G_dn=medfilt1(G_period,'omitnan');

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

H_seg=reshape(H_dn,N_seg,[]);
D_seg=reshape(D_dn,N_seg,[]);
Z_seg=reshape(Z_dn,N_seg,[]);
G_seg=reshape(G_dn,N_seg,[]);
%% %Night time data in LT

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
%% %Diff, PSD and power ratio calculation
Z_diff=diff(Z_dn);
Z_diff(length(time_vec))=NaN;
G_diff=diff(G_dn);
G_diff(length(time_vec))=NaN;

j=1;
for i=1:length(Z_dn)/86400
    Z_filt(j:j+86399,1)=filter1('bp',Z_dn(j:j+86399,1),'fc',[0.01 0.1]);
    G_filt(j:j+86399,1)=filter1('bp',G_dn(j:j+86399,1),'fc',[0.01 0.1]);
    j=j+86400;
end

n_win=1800;
n_ovrlp=0.5*n_win;
n_fft=1800;
fs=1;

bin001=ceil(((n_fft/2+1)/0.5)*0.01);
bin01=ceil(((n_fft/2+1)/0.5)*0.1);

% [H_wlch,f_wlch]=pwelch(H_seg,hamming(n_win),n_ovrlp,n_fft,fs);
% [D_wlch,f_wlch]=pwelch(D_seg,hamming(n_win),n_ovrlp,n_fft,fs);
[Z_wlch,f_wlch]=pwelch(Z_seg,hamming(n_win),n_ovrlp,n_fft,fs);
[G_wlch,f_wlch]=pwelch(G_seg,hamming(n_win),n_ovrlp,n_fft,fs);

% H_ulf=H_wlch(bin001:bin05,:);
% D_ulf=D_wlch(bin001:bin05,:);
Z_ulf=Z_wlch(bin001:bin01,:);
G_ulf=G_wlch(bin001:bin01,:);
f_ulf=f_wlch(bin001:bin01,:);

% H_fgap=fillgaps(H_ulf,20,10);
% D_fgap=fillgaps(D_ulf,20,10);
Z_fgap=fillgaps(Z_ulf,20,10);
G_fgap=fillgaps(G_ulf,20,10);

for i=1:size(Z_ulf,2)

%     if sum(isnan(H_fgap(:,i)))==0
%         fitH=fit(f_ulf,H_fgap(:,i),'exp1');
%         cfitH=coeffvalues(fitH);
%         fittedH=cfitH(1)*exp(cfitH(2)*f_ulf);
%         H_fit(:,i)=H_fgap(:,i)-fittedH;
%     end
    
%     if sum(isnan(D_fgap(:,i)))==0
%         fitD=fit(f_ulf,D_fgap(:,i),'exp1');
%         cfitD=coeffvalues(fitD);
%         fittedD=cfitD(1)*exp(cfitD(2)*f_ulf);
%         D_fit(:,i)=D_fgap(:,i)-fittedD;
%     end

    if sum(isnan(G_fgap(:,i)))==0
        fitG=fit(f_ulf,G_fgap(:,i),'exp1');
        cfitG=coeffvalues(fitG);
        fittedG=cfitG(1)*exp(cfitG(2)*f_ulf);
        G_fit(:,i)=G_fgap(:,i)-fittedG;
    end

    if sum(isnan(Z_fgap(:,i)))==0
        fitZ=fit(f_ulf,Z_fgap(:,1),'exp1');
        cfitZ=coeffvalues(fitZ);
        fittedZ=cfitZ(1)*exp(cfitZ(2)*f_ulf);
        Z_fit(:,i)=Z_fgap(:,i)-fittedZ;
    end
    
end

G_fitmm=movmean(G_fit,5);
Z_fitmm=movmean(Z_fit,5);

bin_1=min((find(f_1==round(f_ulf,3))));
bin_2=max((find(f_2==round(f_ulf,3))));

G_trans=G_fitmm(bin_1:bin_2,:);
Z_trans=Z_fitmm(bin_1:bin_2,:);

[Z_max,Z_i]=max(Z_trans);
Z_norm=(Z_max-min(Z_max))/(max(Z_max)-min(Z_max));

for i=1:length(Z_max)
    G_max(1,i)=(G_trans(Z_i(1,i),i));
end
G_norm=(G_max-min(G_max))/(max(G_max)-min(G_max));

ZG_max=Z_norm./G_norm;

max_f=(((f_2-f_1)/size(Z_trans,1))*Z_i)+f_1;
%% %Acquiring ap and Dst indices
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
else
    ap=AP;
    Dst=DST;
end
%% %Figure 1

title_sup=sprintf('%s station, %s',station,datetime(mean(welchnum_vec),'Format','MMM yyyy'));

f1=figure(1);
stitle=suptitle(sprintf('%s',title_sup));
stitle_pos =get(stitle,'position');
stitle_pos=[stitle_pos(1) stitle_pos(2)+0.01 stitle_pos(3)];
set(stitle,'position',stitle_pos);

spA=subplot(4,2,1)
plot(welchnum_vec,ap,'k');
hold on
plot(welchnum_vec,Dst,'r');
plot(welchnum_vec,ones(length(welchnum_vec),1)*(50),'k--');
plot(welchnum_vec,ones(length(welchnum_vec),1)*(-50),'r--');
hold off
ylabel('Dst               ap');
xlim([min(welchnum_vec) max(welchnum_vec)])
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
        label=sprintf('%s M%0.1f %.0fKM',string(datetime(j,'Format','dd/MM/yyyy')),EQ_sel(i,3),EQ_sel(i,5));
        text(j,EQ_sel(i,4),label,'VerticalAlignment','top','HorizontalAlignment','left')
    end
    hold on
end
hold off
xlim([min(welchnum_vec) max(welchnum_vec)])
ylabel('K_{LS}');
title('Geomagnetic (left) & local seismicity (right) indices');

spB=subplot(4,2,2)
plot(welchnum_vec,ap,'k');
hold on
plot(welchnum_vec,Dst,'r');
plot(welchnum_vec,ones(length(welchnum_vec),1)*(50),'k--');
plot(welchnum_vec,ones(length(welchnum_vec),1)*(-50),'r--');
hold off
ylabel('Dst               ap');
xlim([min(welchnum_vec) max(welchnum_vec)])
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
        label=sprintf('%s M%0.1f %.0fKM',string(datetime(j,'Format','dd/MM/yyyy')),EQ_sel(i,3),EQ_sel(i,5));
        text(j,EQ_sel(i,4),label,'VerticalAlignment','top','HorizontalAlignment','left')
    end
    hold on
end
hold off
xlim([min(welchnum_vec) max(welchnum_vec)])
ylabel('K_{LS}');
title('Geomagnetic (left) & local seismicity (right) indices');
 
spC=subplot(4,2,3)
plot(time_vec,G_dn);
ylabel('G_{raw}')
xlim([min(time_vec) max(time_vec)])
title('Raw G component');

spD=subplot(4,2,4)
plot(time_vec,Z_dn);
ylabel('Z_{raw}')
xlim([min(time_vec) max(time_vec)])
title('Raw Z component');

spE=subplot(4,2,5)
plot(time_vec,G_diff);
ylabel('G_{diff}')
ylim([-3 3])
xlim([min(time_vec) max(time_vec)]);
title('Diff(G)');

spF=subplot(4,2,6)
plot(time_vec,Z_diff);
ylabel('Z_{diff}')
ylim([-3 3])
xlim([min(time_vec) max(time_vec)]);
title('Diff(Z)');

spG=subplot(4,2,7)
plot(time_vec,G_filt);
ylabel('G_{filtered}')
xlim([min(time_vec) max(time_vec)]);
title('G at \Deltaf = 0.01 - 0.1 Hz');

spH=subplot(4,2,8)
plot(time_vec,Z_filt);
ylabel('Z_{filtered}')
xlim([min(time_vec) max(time_vec)]);
title('Z at \Deltaf = 0.01 - 0.1 Hz');

linkaxes([spA,spC,spE,spG,spB,spD,spF,spH],'x')
linkaxes([spG,spH],'y')
set(gcf, 'Position', get(0, 'Screensize'));
%% %Figure 2
title_sup=sprintf('%s station, %s',station,datetime(mean(welchnum_vec),'Format','MMM yyyy'));

f2=figure(2);
stitle=suptitle(sprintf('%s',title_sup));
stitle_pos =get(stitle,'position');
stitle_pos=[stitle_pos(1) stitle_pos(2)+0.01 stitle_pos(3)];
set(stitle,'position',stitle_pos);

sp1=subplot(4,1,1)
plot(welchnum_vec,ap,'k');
hold on
plot(welchnum_vec,Dst,'r');
plot(welchnum_vec,ones(length(welchnum_vec),1)*(50),'k--');
plot(welchnum_vec,ones(length(welchnum_vec),1)*(-50),'r--');
hold off
ylabel('Dst               ap');
xlim([min(welchnum_vec) max(welchnum_vec)])
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
xlim([min(welchnum_vec) max(welchnum_vec)])
ylabel('K_{LS}');
title('Geomagnetic (left) & local seismicity (right) indices');

sp2=subplot(4,1,2)
spectrogram(G_dn,'yaxis',hamming(n_win),n_ovrlp,n_fft,1)
caxis([0 30])
ylim([0 60])
colormap jet
colorbar('off')
colorbar('east')
title('Spectrogram of G');

sp3=subplot(4,1,3)
spectrogram(Z_dn,'yaxis',hamming(n_win),n_ovrlp,n_fft,1)
caxis([0 30])
ylim([0 60])
colormap jet
colorbar('off')
colorbar('east')
title('Spectrogram of Z');

sp4=subplot(4,1,4)
spln=spline(1:length(welchnum_vec),max_f,1:1:length(welchnum_vec));
plot(welchnum_vec,max_f,'ko','MarkerSize',2,'MarkerFaceColor','k')
hold on
plot(welchnum_vec(1:1:end),spln,'r');
hold off
ylabel('Frequency (Hz)')
xlim([min(welchnum_vec) max(welchnum_vec)])
ylim([0.8*min(max_f) 1.2*max(max_f)])
spln=spline(1:length(welchnum_vec),max_f,1:1:length(welchnum_vec));
for i=0:size(welchnum_vec,1)/(86400/N_seg)
    xpatch=[t_start/24+i t_end/24+i t_end/24+i t_start/24+i];
    ypatch=[-1e+5 -1e+5 1e+5 1e+5];
    patch(xpatch,ypatch,'k','FaceAlpha',0.05,'EdgeAlpha',0);
end
yyaxis right
plot(welchnum_vec,ZG_max,'LineWidth',0.8)
ylabel('Z/G_{Max power}')
xlim([min(welchnum_vec) max(welchnum_vec)])
ylim([0.5*min(ZG_max) 1.2*max(ZG_max)])
title(sprintf('Dominant frequencies (left) & Maximum power of Z/H (right) in %0.3f - %.3f Hz',f_1,f_2));
[val ind]=sort(ZG_max,'descend');
maxval=val(1:3);
maxind=ind(1:3);
for i=1:3
    label=sprintf('%0.3f Hz',max_f(1,maxind(1,i)));
    text(welchnum_vec(maxind(1,i)),maxval(1,i),label,'VerticalAlignment','bottom','HorizontalAlignment','left')
end

linkaxes([sp1,sp4],'x')
set(gcf, 'Position', get(0, 'Screensize'));
%% 
% %% %Figure 3
% bpf3=input('Narrower frequency range? f1 ')
% bpf4=input('Narrower frequency range? f2 ')
% 
% H_filt=filter1('bp',H_dn,'fc',[bpf3 bpf4]);
% Z_filt=filter1('bp',Z_dn,'fc',[bpf3 bpf4]);
% 
% f3=figure(3);
% stitle=suptitle(sprintf('%s',title_sup));
% stitle_pos =get(stitle,'position');
% stitle_pos=[stitle_pos(1) stitle_pos(2)+0.01 stitle_pos(3)];
% set(stitle,'position',stitle_pos);
% 
% subplot(3,2,[1 2])
% plot(welchnum_vec,ap,'k');
% hold on
% plot(welchnum_vec,Dst,'r');
% plot(welchnum_vec,ones(length(welchnum_vec),1)*(50),'k--');
% plot(welchnum_vec,ones(length(welchnum_vec),1)*(-50),'r--');
% hold off
% ylabel('Dst               ap');
% xlim([min(welchnum_vec) max(welchnum_vec)])
% if max(ap)>50 || min(Dst)<-50
%     ylim([min(Dst)-10 max(ap)+10])
% else
%     ylim([-60 60])
% end
% yyaxis right
% set(gca, 'YColor', 'b')
% for i=1:size(EQ_sel,1)
%     j=datetime(datevec(EQ_sel(i,1)),'Format','dd/MM/yyyy HH:mm:ss');
%     plot(j,EQ_sel(i,4),'bo','MarkerFaceColor','b','MarkerSize',5);
%     if EQ_sel(i,4)==max(EQ_sel(:,4))
%         label=sprintf('%s M%0.1f %.0fKM',string(datetime(j,'Format','dd/MM/yyyy')),EQ_sel(i,3),EQ_sel(i,5));
%         text(j,EQ_sel(i,4),label,'VerticalAlignment','top','HorizontalAlignment','left')
%     end
%     hold on
% end
% hold off
% xlim([min(welchnum_vec) max(welchnum_vec)])
% ylabel('K_{LS}');
% title('Geomagnetic (left) & local seismicity (right) indices');
% 
% subplot(3,2,3)
% plot(time_vec,H_filt,'r')
% ylabel('H_{filtered}');
% xlim([min(time_vec) max(time_vec)])
% title(sprintf('Bandpass filtered of H (%0.3f - %.3f Hz)',bpf3,bpf4));
% 
% subplot(3,2,4)
% plot(time_vec,Z_filt,'r')
% ylabel('Z_{filtered}');
% xlim([min(time_vec) max(time_vec)])
% title(sprintf('Bandpass filtered of Z (%0.3f - %.3f Hz)',bpf3,bpf4));
% 
% ZH=Z_dn./H_dn;
% ZH_filt=filter1('bp',ZH,'fc',[bpf3 bpf4]);
% subplot(3,2,[5 6])
% plot(time_vec,ZH_filt)
% ylabel('Z/H_{filtered}');
% xlim([min(time_vec) max(time_vec)])
% title(sprintf('Bandpass filtered of Z/H (%0.3f - %.3f Hz)',bpf3,bpf4));
% 
% set(gcf, 'Position', get(0, 'Screensize'));
%% %Saving figures
MMM_yyyy=upper(string(datetime(mean(welchnum_vec),'Format','MMMyyyy')));

today=char(datetime('today','Format','dd-MM-yyyy'));
today=strcat(today,'\BMKG_5\');
path1='D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\';
mkdir(fullfile(path1,today));

path2=strcat('D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\',today);

% figname1=char(strcat(stn,'-',MMM_yyyy,' (RAW)',addtitle));
% saveas(f1,fullfile(path2,figname1),'tiff');
% saveas(f1,fullfile(path2,figname1),'fig');

fq_1=regexprep(string(f_1),'[.]','');
fq_2=regexprep(string(f_2),'[.]','');
figname2=char(strcat(stn,'-',MMM_yyyy,'-',fq_1,'-',fq_2,' (MAXPOW)',addtitle));
saveas(f2,fullfile(path2,figname2),'tiff');

% fq_3=regexprep(string(bpf3),'[.]','');
% fq_4=regexprep(string(bpf4),'[.]','');
% figname3=char(strcat(stn,'-',MMM_yyyy,'-',fq_3,'-',fq_4,' (BANDPASS)',addtitle));
% saveas(f3,fullfile(path2,figname3),'tiff');
% saveas(f3,fullfile(path2,figname3),'fig');
