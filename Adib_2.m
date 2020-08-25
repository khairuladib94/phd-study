%Input---------------------------------------------------------------------
stn='KTB';                           %Abbreviation of station name
station='Kototabang';                %Station full name
mag_min=5.0;                         %Minimum magnitude of earthquakes to be considered
dis_max=550;                         %Maximum epicentral distance from the station
range='FULL';                         %Write 'FULL' if desired period is the whole year. Otherwise, leave blank

day_start=1;                         %Insert start date and end date
month_start=6;
year_start=2010;

day_end=31;
month_end=12;
year_end=2010;

year=year_start;

%--------------------------------------------------------------------------

DOY_start=day(datetime([year_start month_start day_start]),'dayofyear');
DOY_end=day(datetime([year_end month_end day_end]),'dayofyear');

days_num=DOY_end-DOY_start+1;
days=DOY_start:1:DOY_end;

if range=='FULL'
    year=min(datevec(UT1m));
    year=year(1,1);
    DOY_start=1;
    if rem(year,4)==0
        DOY_end=366;
    end
    if rem(year,4)~=0
        DOY_end=365;
    end
    days_num=DOY_end;
    days=DOY_start:1:DOY_end;
end

%Setting station number
if stn=='TGG'
    stn_num=1;
end
if stn=='MUT'
    stn_num=2;
end
if stn=='LGZ'
    stn_num=3;
end
if stn=='CEB'
    stn_num=4;
end
if stn=='CDO'
    stn_num=5;
end
if stn=='DAV'
    stn_num=6;
end
if stn=='GSI'
    stn_num=7;
end
if stn=='SCN'
    stn_num=8;
end
if stn=='LWA'
    stn_num=9;
end
if stn=='PTN'
    stn_num=10;
end
if stn=='MND'
    stn_num=11;
end
if stn=='BIK'
    stn_num=12;
end
if stn=='JYP'
    stn_num=13;
end
if stn=='PRP'
    stn_num=14;
end
if stn=='KPG'
    stn_num=15;
end
if stn=='KTB'
    stn_num=16;
end
    
%Getting station coordinate
stn_latlon=[stn_IndPhi(stn_num,2) stn_IndPhi(stn_num,3)];

%Earthquake----------------------------------------------------------------

%Building selected earthquakes table
EQ_sel=NaN(100,5);
i=1;
for j=1:length(EQ_IndPhi1117)
    if EQ_IndPhi1117(j,1)==year
        if  EQ_IndPhi1117(j,10)>=mag_min
            EQ_latlon=[EQ_IndPhi1117(j,7) EQ_IndPhi1117(j,8)];
            EQ_dis=deg2km(distance('gc',stn_latlon,EQ_latlon));
            EQ_date=datenum([EQ_IndPhi1117(j,1) EQ_IndPhi1117(j,2) EQ_IndPhi1117(j,3)]);
            EQ_mag=EQ_IndPhi1117(j,10);
            EQ_DOY=day(datetime([EQ_IndPhi1117(j,1) EQ_IndPhi1117(j,2) EQ_IndPhi1117(j,3)]),'dayofyear');
            EQ_Ks=((1+EQ_dis^(-EQ_mag/2))^(-2.33))+((10^(0.75*EQ_mag))/(10*EQ_dis));  
 
            if EQ_dis<=dis_max
                EQ_sel(i,1)=EQ_date;
                EQ_sel(i,2)=EQ_DOY;
                EQ_sel(i,3)=EQ_dis;
                EQ_sel(i,4)=EQ_mag;
                EQ_sel(i,5)=EQ_Ks;
                i=i+1;
            end
        end
    end
end

EQ_sel(any(isnan(EQ_sel),2),:)=[];
if size(EQ_sel,1)>1
    EQ_sel=flip(EQ_sel);
end

%Plotting seismic index chart

figure(1)

subplot(5,1,1)
plot(EQ_sel(:,2),EQ_sel(:,5),'rx')
xlim([DOY_start DOY_end]);
ylim([0 max(EQ_sel(:,5))+10])
xlabel('Day of year');
ylabel('Seismic Index (K_{s})');
title(sprintf('Polarisation ratio | %s station, %d',station,year));

display('Done p.1')

%Geomagnetic data----------------------------------------------------------

H_period=H(86400*DOY_start-86399:86400*DOY_end);
D_period=D(86400*DOY_start-86399:86400*DOY_end);
Z_period=Z(86400*DOY_start-86399:86400*DOY_end);

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

display('Done p.2')

%Reshaping into daily data
H_day=reshape(H_dn,86400,[]);
D_day=reshape(D_dn,86400,[]);
Z_day=reshape(Z_dn,86400,[]);

%Nighttime data 2100 - 0200 hrs LT. Insert time in UT
t_start=13;      %Phillipines/Indonesia timezone -> UT+8
t_stop=18;

t_start=t_start*3600;      
t_stop=t_stop*3600;
t_int=t_stop-t_start;

H_night=H_day(t_start:t_stop-1,:);
D_night=D_day(t_start:t_stop-1,:);
Z_night=Z_day(t_start:t_stop-1,:);

%Remerging data
H_merge=reshape(H_night,[],1);
D_merge=reshape(D_night,[],1);
Z_merge=reshape(Z_night,[],1);

%DOY ticks setup
DOY_num=t_int:t_int:length(H_merge);

%Plotting geomagnetic raw nighttime data
figure(2)

subplot(3,1,1);

plot(H_merge);
set(gca,'XTick',[DOY_num],'XTickLabel',{days});
xlabel('Day of year');
ylabel('H_{raw}');
xlim([1 length(H_merge)]);
title(sprintf('HDZ nighttime raw data | %s station, %d',station,year));

subplot(3,1,2);
plot(D_merge);
set(gca,'XTick',[DOY_num],'XTickLabel',{days});
xlabel('Day of year');
ylabel('D_{raw}');
xlim([1 length(H_merge)]);

subplot(3,1,3);
plot(Z_merge);
set(gca,'XTick',[DOY_num],'XTickLabel',{days});
xlabel('Day of year');
ylabel('Z_{raw}');
xlim([1 length(H_merge)]);


%Identifying segments of data by excluding gaps that span for 1 day or longer
H_nan=NaN(numel(H_merge),3);
H_nan(:,1)=1:numel(H_merge);
H_nan(:,2)=isnan(H_merge);

H_seg=NaN(100,1);

j=1;
k=1;
for i=1:numel(H_merge)
    if i~=1
        if H_nan(i,2)==1 && H_nan(i-1,2)==1
            k=k+1;
            if H_nan(i+1,2)==0
                if k>=9000
                    H_nan(i,3)=k;
                    H_seg(j,1)=H_nan(i,1)-k+1;
                    j=j+1;
                    H_seg(j,1)=H_nan(i,1);
                    j=j+1;
                    k=1;
                end
            end
        end
    end
end
H_seg(any(isnan(H_seg),2),:)=[];


%Setting up bandpass filter parameters
f_s=1;
filtN=6;
filtfL=0.01;
filtfH=0.1;

[b,a]=butter(filtN,[filtfL filtfH]/(f_s/2),'bandpass');

%Applying bandpass filter and gap filling on each segment
H_fg=NaN(numel(H_merge),1);
H_filt=NaN(numel(H_merge),1);

D_fg=NaN(numel(D_merge),1);
D_filt=NaN(numel(D_merge),1);

Z_fg=NaN(numel(Z_merge),1);
Z_filt=NaN(numel(Z_merge),1);

j=1;
for i=1:numel(H_seg)
    
    if j>numel(H_seg)
        break;
    end
    
    if j==1
        H_fg(1:H_seg(j)-1)=fillgaps(H_merge(1:H_seg(j)-1));
        H_filt(1:H_seg(j)-1)=filter(b,a,H_fg(1:H_seg(j)-1));
        
        D_fg(1:H_seg(j)-1)=fillgaps(D_merge(1:H_seg(j)-1));
        D_filt(1:H_seg(j)-1)=filter(b,a,D_fg(1:H_seg(j)-1));
        
        Z_fg(1:H_seg(j)-1)=fillgaps(Z_merge(1:H_seg(j)-1));
        Z_filt(1:H_seg(j)-1)=filter(b,a,Z_fg(1:H_seg(j)-1));
    end
    
    if j~=1 && j~=numel(H_seg)
        H_fg(H_seg(j)+1:H_seg(j+1)-1)=fillgaps(H_merge(H_seg(j)+1:H_seg(j+1)-1));
        H_filt(H_seg(j)+1:H_seg(j+1)-1)=filter(b,a,H_fg(H_seg(j)+1:H_seg(j+1)-1));
        
        D_fg(H_seg(j)+1:H_seg(j+1)-1)=fillgaps(D_merge(H_seg(j)+1:H_seg(j+1)-1));
        D_filt(H_seg(j)+1:H_seg(j+1)-1)=filter(b,a,D_fg(H_seg(j)+1:H_seg(j+1)-1));
        
        Z_fg(H_seg(j)+1:H_seg(j+1)-1)=fillgaps(Z_merge(H_seg(j)+1:H_seg(j+1)-1));
        Z_filt(H_seg(j)+1:H_seg(j+1)-1)=filter(b,a,Z_fg(H_seg(j)+1:H_seg(j+1)-1));
        
        j=j+1;
    end
    
    if j==numel(H_seg)
        H_fg(H_seg(j)+1:numel(H_merge))=fillgaps(H_merge(H_seg(j)+1:numel(H_merge)));
        H_filt(H_seg(j)+1:numel(H_merge))=filter(b,a,H_fg(H_seg(j)+1:numel(H_merge)));
        
        D_fg(H_seg(j)+1:numel(H_merge))=fillgaps(D_merge(H_seg(j)+1:numel(H_merge)));
        D_filt(H_seg(j)+1:numel(H_merge))=filter(b,a,D_fg(H_seg(j)+1:numel(H_merge)));
       
        Z_fg(H_seg(j)+1:numel(H_merge))=fillgaps(Z_merge(H_seg(j)+1:numel(H_merge)));
        Z_filt(H_seg(j)+1:numel(H_merge))=filter(b,a,Z_fg(H_seg(j)+1:numel(H_merge)));
    end
    
    j=j+1;
    
end

display('Done p.3')

% %Plotting geomagnetic gap-filled, filtered data
figure(3)

subplot(3,1,1);
plot(H_filt);
set(gca,'XTick',[DOY_num],'XTickLabel',{days});
xlabel('Day of year');
ylabel('H_{filtered}');
xlim([1 length(H_filt)]);
title(sprintf('HDZ nighttime gap-filled, bandpass-filtered data | %s station, %d',station,year));

subplot(3,1,2);
plot(D_filt);
set(gca,'XTick',[DOY_num],'XTickLabel',{days});
xlabel('Day of year');
ylabel('D_{filtered}');
xlim([1 length(H_filt)]);

subplot(3,1,3);
plot(Z_filt);
set(gca,'XTick',[DOY_num],'XTickLabel',{days});
xlabel('Day of year');
ylabel('Z_{filtered}');
xlim([1 length(H_filt)]);

%Reshaping into daily data
H_day2=reshape(H_filt,t_int,[]);
D_day2=reshape(D_filt,t_int,[]);
Z_day2=reshape(Z_filt,t_int,[]);

%Trimming out the first one hour per day data to remove initial transient
H_trim=H_day2(3601:t_int,:);    
D_trim=D_day2(3601:t_int,:);
Z_trim=Z_day2(3601:t_int,:);

%Remerging the data
H_merge2=reshape(H_trim,[],1);
D_merge2=reshape(D_trim,[],1);
Z_merge2=reshape(Z_trim,[],1);

%Calculating mean and standard deviation for H field
H_sig=std(H_merge2,'omitnan');
H_mu=mean(H_merge2,'omitnan');

H_mupsig=H_mu+3*H_sig;
H_mumsig=H_mu-3*H_sig;

t_int2=t_int-3600;

% %Removing any leftover initial transient
% for i=1:length(H_merge2)
%     
%     if H_merge2(i)>H_mupsig||H_merge2(i)<H_mumsig
%         i_rem=i/t_int2;
%         i_rem=floor(i_rem);
%         i_start=i_rem*t_int2+1;
%         i_stop=i_start+t_int2-1;
%         
%         H_merge2(i_start:i_stop)=NaN;
%         D_merge2(i_start:i_stop)=NaN;
%         Z_merge2(i_start:i_stop)=NaN;
%     end
%     
% end

display('Done p.4')

DOY_num2=t_int2:t_int2:length(H_merge2);

%Plotting geomagnetic initial transient trimmed data
figure(4)

subplot(3,1,1);
plot(H_merge2);
set(gca,'XTick',[DOY_num2],'XTickLabel',{days});
xlabel('Day of year');
ylabel('H_{filtered}');
xlim([1 length(H_merge2)]);
title(sprintf('HDZ nighttime initial-transient trimmed | %s station, %d',station,year));

subplot(3,1,2);
plot(D_merge2);
set(gca,'XTick',[DOY_num2],'XTickLabel',{days});
xlabel('Day of year');
ylabel('D_{filtered}');
xlim([1 length(H_merge2)]);

subplot(3,1,3);
plot(Z_merge2);
set(gca,'XTick',[DOY_num2],'XTickLabel',{days});
xlabel('Day of year');
ylabel('Z_{filtered}');
xlim([1 length(H_merge2)]);

%Reshaping into daily data
H_day3=reshape(H_merge2,t_int2,[]);
D_day3=reshape(D_merge2,t_int2,[]);
Z_day3=reshape(Z_merge2,t_int2,[]);

%Calculating daily mean
H_daymu=mean(H_day3,1,'omitnan');
D_daymu=mean(D_day3,1,'omitnan');
Z_daymu=mean(Z_day3,1,'omitnan');

%Remerging the daily mean data
H_merge3=reshape(H_daymu,[],1);
D_merge3=reshape(D_daymu,[],1);
Z_merge3=reshape(Z_daymu,[],1);

%Plotting geomagnetic daily mean data
figure(5)

subplot(3,1,1);
plot(days,H_merge3);
xlabel('Day of year');
ylabel('H_{filtered}');
xlim([DOY_start DOY_end]);

subplot(3,1,2);
plot(days,D_merge3);
xlabel('Day of year');
ylabel('D_{filtered}');
xlim([DOY_start DOY_end]);

subplot(3,1,3);
plot(days,Z_merge3);
xlabel('Day of year');
ylabel('Z_{filtered}');
xlim([DOY_start DOY_end]);

%Mean for the whole period
H_allmu=mean(H_merge2,1,'omitnan');
D_allmu=mean(D_merge2,1,'omitnan');
Z_allmu=mean(Z_merge2,1,'omitnan');

%Standard deviation for whole period
H_allsig=std(H_merge2,0,1,'omitnan');
D_allsig=std(D_merge2,0,1,'omitnan');
Z_allsig=std(Z_merge2,0,1,'omitnan');

%Normalisation
H_norm=(H_daymu-H_allmu)/H_allsig;
D_norm=(D_daymu-D_allmu)/D_allsig;
Z_norm=(Z_daymu-Z_allmu)/Z_allsig;

%H and D exceeding 0.1 only and polarisation ratio
ZH2=NaN(1,days_num);
ZD2=NaN(1,days_num);

for i=1:days_num
    if H_norm(1,i)>0.1
        ZH2(1,i)=Z_merge3(i,1)/H_merge3(i,1);
    else
        ZH2(1,i)=0;
    end
    if D_norm(1,i)>0.1
        ZD2(1,i)=Z_merge3(i,1)/D_merge3(i,1);
    else
        ZD2(1,i)=0;
    end
end

display('Done p.5')

%Calculating polarisation ratio
ZH=Z_merge3./H_merge3;
ZD=Z_merge3./D_merge3;
ZG=Z_merge3./(sqrt(H_merge3.^2+D_merge3.^2));

ZH_mm=movmean(ZH,5,'omitnan');
ZD_mm=movmean(ZD,5,'omitnan');
ZG_mm=movmean(ZG,5,'omitnan');

%Plotting polarisation ratio
figure(1)

subplot(5,1,2);

plot(days,ZH,days,ZH_mm,'r');
xlabel('Day of year');
ylabel('Z/H');
xlim([DOY_start DOY_end]);

subplot(5,1,3);
plot(days,ZD,days,ZD_mm,'r');
xlabel('Day of year');
ylabel('Z/D');
xlim([DOY_start DOY_end]);

subplot(5,1,4);
plot(days,ZG,days,ZG_mm,'r');
xlabel('Day of year');
ylabel('Z/G');
xlim([DOY_start DOY_end]);

%Getting Kp and Dst indices
Kp=zeros(days_num,1);
Dst=zeros(days_num,1);
j=1;
for i=1:length(index_1117)
    if index_1117(i,1)==year && index_1117(i,2)>=DOY_start
        Kp(j,1)=index_1117(i,4);
        Dst(j,1)=index_1117(i,5);
        j=j+1;
        
        if j==days_num+1
            break;
        end
    end
end

%Plotting Kp and Dst indices
figure(1)

subplot(5,1,5)
bar(days,Kp,'blue');
hold on
bar(days,Dst,'red');
plot(days,ones(days_num,1)*(32),'b--');
plot(days,ones(days_num,1)*(-50),'r--');
hold off
xlabel('Day of year');
ylabel('Dst     \SigmaK_{p}');
legend('\SigmaK_{p}','Dst','location','southeast');
xlim([DOY_start DOY_end]);
ylim([-60 40]);

display('Done p.6')
