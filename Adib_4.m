%Input---------------------------------------------------------------------
stn='KTB';                           %Abbreviation of station name
station='Kotatabang';                %Station full name
mag_min=5.0;                         %Minimum magnitude of earthquakes to be considered
dis_max=300;                         %Maximum epicentral distance from the station
range='FULL';                         %Write 'FULL' if desired period is the whole year. Otherwise, leave blank
nfft=2048;

day_start=1;                         %Insert start date and end date
month_start=9;
year_start=2010;

day_end=24;
month_end=11;
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
    else 
        DOY_end=365;
    end
    days_num=DOY_end;
    days=DOY_start:1:DOY_end;
end

%Setting station number
if stn=='TGG'
    stn_num=1;
elseif stn=='MUT'
    stn_num=2;
elseif stn=='LGZ'
    stn_num=3;
elseif stn=='CEB'
    stn_num=4;
elseif stn=='CDO'
    stn_num=5;
elseif stn=='DAV'
    stn_num=6;
elseif stn=='GSI'
    stn_num=7;
elseif stn=='SCN'
    stn_num=8;
elseif stn=='LWA'
    stn_num=9;
elseif stn=='PTN'
    stn_num=10;
elseif stn=='MND'
    stn_num=11;
elseif stn=='BIK'
    stn_num=12;
elseif stn=='JYP'
    stn_num=13;
elseif stn=='PRP'
    stn_num=14;
elseif stn=='KPG'
    stn_num=15;
elseif stn=='KTB'
    stn_num=16;
end
    
%Getting station coordinate
stn_latlon=[stn_IndPhi(stn_num,2) stn_IndPhi(stn_num,3)];

%Earthquake----------------------------------------------------------------

%Building selected earthquakes table
EQ_sel=NaN(100,7);
i=1;
for j=1:length(EQ_IndPhi1117)
    if EQ_IndPhi1117(j,1)==year
        if  EQ_IndPhi1117(j,10)>=mag_min
            EQ_latlon=[EQ_IndPhi1117(j,7) EQ_IndPhi1117(j,8)];
            EQ_dis=deg2km(distance('gc',stn_latlon,EQ_latlon));
            EQ_date=datenum([EQ_IndPhi1117(j,1) EQ_IndPhi1117(j,2) EQ_IndPhi1117(j,3)]);
            EQ_mag=EQ_IndPhi1117(j,10);
            EQ_DOY=day(datetime([EQ_IndPhi1117(j,1) EQ_IndPhi1117(j,2) EQ_IndPhi1117(j,3)]),'dayofyear');
            EQ_Ks=((10^(0.75*EQ_mag))/(10*EQ_dis))*((1+EQ_dis*10^(-EQ_mag/2))^(-2.33));  
            EQ_depth=EQ_IndPhi1117(j,9);
            EQ_f=(2*600)/(((1000*EQ_depth)^2)*1.2566e-6*2*pi);
            
            if EQ_dis<=dis_max
                EQ_sel(i,1)=EQ_date;
                EQ_sel(i,2)=EQ_DOY;
                EQ_sel(i,3)=EQ_dis;
                EQ_sel(i,4)=EQ_mag;
                EQ_sel(i,5)=EQ_Ks;
                EQ_sel(i,6)=EQ_depth;
                EQ_sel(i,7)=EQ_f;
                i=i+1;
            end
        end
    end
end

EQ_sel(any(isnan(EQ_sel),2),:)=[];
if size(EQ_sel,1)>1
    EQ_sel=flip(EQ_sel);
end

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

%Reshaping into daily data
H_day=reshape(H_dn,86400,[]);
D_day=reshape(D_dn,86400,[]);
Z_day=reshape(Z_dn,86400,[]);

%Nighttime data 2200 - 0200 hrs LT. Insert time in UT
t_start=14;      %Phillipines/Indonesia timezone -> UT+8
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

%Splitting data into 30mins intervals. 1 day=8 segments
H_30min=reshape(H_merge,1800,[]);
D_30min=reshape(D_merge,1800,[]);
Z_30min=reshape(Z_merge,1800,[]);

f_A1=0.008;
f_A2=0.012;
f_B1=0.08;
f_B2=0.12;
f_C1=0.45;
f_C2=0.49;

% f_A1=0.01;
% f_A2=0.1;
% f_B1=0.01;
% f_B2=0.015;
% f_C1=0.02;
% f_C2=0.05;
% 
% f_A1=0.03;
% f_A2=0.04;
% f_B1=0.04;
% f_B2=0.05;
% f_C1=0.05;
% f_C2=0.06;

% f_A1=0.0005;
% f_A2=0.0010;
% f_B1=0.0010;
% f_B2=0.0015;
% f_C1=0.0015;
% f_C2=0.0020;

%---------------------------------------------------------

H_wmu=NaN(8,days_num);
D_wmu=NaN(8,days_num);
Z_wmu=NaN(8,days_num);

bin_1=round(((nfft/2+1)/0.5)*f_A1);
bin_2=round(((nfft/2+1)/0.5)*f_A2);

f_s=1;
j=1;
k=1;

for i=1:size(H_30min,2)
        [H_welch,f]=pwelch(H_30min(:,i),1800,[],nfft,f_s);
        H_welch=abs(H_welch(bin_1:bin_2));
        H_wmu(j,k)=mean(H_welch);
        
        [D_welch,f]=pwelch(D_30min(:,i),1800,[],nfft,f_s);
        D_welch=abs(D_welch(bin_1:bin_2));
        D_wmu(j,k)=mean(D_welch);
        
        [Z_welch,f]=pwelch(Z_30min(:,i),1800,[],nfft,f_s);
        Z_welch=abs(Z_welch(bin_1:bin_2));
        Z_wmu(j,k)=mean(Z_welch);
        
        j=j+1;
  
        if j>8
            k=k+1;
            j=1;
        end
end

H_daymu=10*log10(mean(H_wmu,1));
D_daymu=10*log10(mean(D_wmu,1));
Z_daymu=10*log10(mean(Z_wmu,1));

%---------------------------------

H_monmu=NaN(1,12);
D_monmu=NaN(1,12);
Z_monmu=NaN(1,12);

H_monsig=NaN(1,12);
D_monsig=NaN(1,12);
Z_monsig=NaN(1,12);

if days_num==365
    days_mon=[31,28,31,30,31,30,31,31,30,31,30,31];
else
    days_mon=[31,29,31,30,31,30,31,31,30,31,30,31];
end

j=1;
k=0;

for i=1:12
    
    k=k+days_mon(i);
    
    H_monmu(1,i)=mean(H_daymu(j:k),'omitnan');
    D_monmu(1,i)=mean(D_daymu(j:k),'omitnan');
    Z_monmu(1,i)=mean(Z_daymu(j:k),'omitnan');
    
    H_monsig(1,i)=std(H_daymu(j:k),'omitnan');
    D_monsig(1,i)=std(D_daymu(j:k),'omitnan');
    Z_monsig(1,i)=std(Z_daymu(j:k),'omitnan');
    
    j=j+days_mon(i);
    
end


%Normalisation

i=1;
j=1;
for k=1:days_num
        
    H_norm(1,k)=(H_daymu(1,k)-H_monmu(1,j))/H_monsig(1,j);
    D_norm(1,k)=(D_daymu(1,k)-D_monmu(1,j))/D_monsig(1,j);
    Z_norm(1,k)=(Z_daymu(1,k)-Z_monmu(1,j))/Z_monsig(1,j);
    
    i=i+1;
    
    if i>days_mon(j)
        j=j+1;
        i=1;
    end
  
end

%----------------------------------

ZG1=Z_norm./(sqrt(H_norm.^2+D_norm.^2));

%mu+/-2sigma
ZG1_mu=mean(ZG1,'omitnan');
ZG1_sig=std(ZG1,'omitnan');
ZG1_mup2sig=ones(1,days_num)*(ZG1_mu+2*ZG1_sig);
ZG1_mum2sig=ones(1,days_num)*(ZG1_mu-2*ZG1_sig);
ZG1_mm=movmean(ZG1,5,'omitnan');

%------------------------------------

H_wmu=NaN(8,days_num);
D_wmu=NaN(8,days_num);
Z_wmu=NaN(8,days_num);

bin_1=round(((nfft/2+1)/0.5)*f_B1);
bin_2=round(((nfft/2+1)/0.5)*f_B2);

f_s=1;
j=1;
k=1;

for i=1:size(H_30min,2)
        [H_welch,f]=pwelch(H_30min(:,i),1800,[],nfft,f_s);
        H_welch=abs(H_welch(bin_1:bin_2));
        H_wmu(j,k)=mean(H_welch);
        
        [D_welch,f]=pwelch(D_30min(:,i),1800,[],nfft,f_s);
        D_welch=abs(D_welch(bin_1:bin_2));
        D_wmu(j,k)=mean(D_welch);
        
        [Z_welch,f]=pwelch(Z_30min(:,i),1800,[],nfft,f_s);
        Z_welch=abs(Z_welch(bin_1:bin_2));
        Z_wmu(j,k)=mean(Z_welch);
        
        j=j+1;
  
        if j>8
            k=k+1;
            j=1;
        end
end

H_daymu=10*log10(mean(H_wmu,1));
D_daymu=10*log10(mean(D_wmu,1));
Z_daymu=10*log10(mean(Z_wmu,1));

%---------------------------------

H_monmu=NaN(1,12);
D_monmu=NaN(1,12);
Z_monmu=NaN(1,12);

H_monsig=NaN(1,12);
D_monsig=NaN(1,12);
Z_monsig=NaN(1,12);

if days_num==365
    days_mon=[31,28,31,30,31,30,31,31,30,31,30,31];
else
    days_mon=[31,29,31,30,31,30,31,31,30,31,30,31];
end

j=1;
k=0;

for i=1:12
    
    k=k+days_mon(i);
    
    H_monmu(1,i)=mean(H_daymu(j:k),'omitnan');
    D_monmu(1,i)=mean(D_daymu(j:k),'omitnan');
    Z_monmu(1,i)=mean(Z_daymu(j:k),'omitnan');
    
    H_monsig(1,i)=std(H_daymu(j:k),'omitnan');
    D_monsig(1,i)=std(D_daymu(j:k),'omitnan');
    Z_monsig(1,i)=std(Z_daymu(j:k),'omitnan');
    
    j=j+days_mon(i);
    
end

%Normalisation

i=1;
j=1;
for k=1:days_num
        
    H_norm(1,k)=(H_daymu(1,k)-H_monmu(1,j))/H_monsig(1,j);
    D_norm(1,k)=(D_daymu(1,k)-D_monmu(1,j))/D_monsig(1,j);
    Z_norm(1,k)=(Z_daymu(1,k)-Z_monmu(1,j))/Z_monsig(1,j);
    
    i=i+1;
    
    if i>days_mon(j)
        j=j+1;
        i=1;
    end
  
end

ZG2=Z_norm./(sqrt(H_norm.^2+D_norm.^2));

%mu+/-2sigma
ZG2_mu=mean(ZG2,'omitnan');
ZG2_sig=std(ZG2,'omitnan');
ZG2_mup2sig=ones(1,days_num)*(ZG2_mu+2*ZG2_sig);
ZG2_mum2sig=ones(1,days_num)*(ZG2_mu-2*ZG2_sig);
ZG2_mm=movmean(ZG2,5,'omitnan');


%----------------------------------------


H_wmu=NaN(8,days_num);
D_wmu=NaN(8,days_num);
Z_wmu=NaN(8,days_num);

bin_1=round(((nfft/2+1)/0.5)*f_C1);
bin_2=round(((nfft/2+1)/0.5)*f_C2);

f_s=1;
j=1;
k=1;

for i=1:size(H_30min,2)
        [H_welch,f]=pwelch(H_30min(:,i),1800,[],nfft,f_s);
        H_welch=abs(H_welch(bin_1:bin_2));
        H_wmu(j,k)=mean(H_welch);
        
        [D_welch,f]=pwelch(D_30min(:,i),1800,[],nfft,f_s);
        D_welch=abs(D_welch(bin_1:bin_2));
        D_wmu(j,k)=mean(D_welch);
        
        [Z_welch,f]=pwelch(Z_30min(:,i),1800,[],nfft,f_s);
        Z_welch=abs(Z_welch(bin_1:bin_2));
        Z_wmu(j,k)=mean(Z_welch);
        
        j=j+1;
  
        if j>8
            k=k+1;
            j=1;
        end
end

H_daymu=10*log10(mean(H_wmu,1));
D_daymu=10*log10(mean(D_wmu,1));
Z_daymu=10*log10(mean(Z_wmu,1));

%---------------------------------

H_monmu=NaN(1,12);
D_monmu=NaN(1,12);
Z_monmu=NaN(1,12);

H_monsig=NaN(1,12);
D_monsig=NaN(1,12);
Z_monsig=NaN(1,12);

if days_num==365
    days_mon=[31,28,31,30,31,30,31,31,30,31,30,31];
else
    days_mon=[31,29,31,30,31,30,31,31,30,31,30,31];
end

j=1;
k=0;

for i=1:12
    
    k=k+days_mon(i);
    
    H_monmu(1,i)=mean(H_daymu(j:k),'omitnan');
    D_monmu(1,i)=mean(D_daymu(j:k),'omitnan');
    Z_monmu(1,i)=mean(Z_daymu(j:k),'omitnan');
    
    H_monsig(1,i)=std(H_daymu(j:k),'omitnan');
    D_monsig(1,i)=std(D_daymu(j:k),'omitnan');
    Z_monsig(1,i)=std(Z_daymu(j:k),'omitnan');
    
    j=j+days_mon(i);
    
end

%Normalisation

i=1;
j=1;
for k=1:days_num
        
    H_norm(1,k)=(H_daymu(1,k)-H_monmu(1,j))/H_monsig(1,j);
    D_norm(1,k)=(D_daymu(1,k)-D_monmu(1,j))/D_monsig(1,j);
    Z_norm(1,k)=(Z_daymu(1,k)-Z_monmu(1,j))/Z_monsig(1,j);
    
    i=i+1;
    
    if i>days_mon(j)
        j=j+1;
        i=1;
    end
  
end

ZG3=Z_norm./(sqrt(H_norm.^2+D_norm.^2));

%mu+/-2sigma
ZG3_mu=mean(ZG3,'omitnan');
ZG3_sig=std(ZG3,'omitnan');
ZG3_mup2sig=ones(1,days_num)*(ZG3_mu+2*ZG3_sig);
ZG3_mum2sig=ones(1,days_num)*(ZG3_mu-2*ZG3_sig);
ZG3_mm=movmean(ZG3,5,'omitnan');


%----------------------------------------------

%Plotting seismic index chart

subplot(5,1,1)
for i=1:size(EQ_sel,1)
    if EQ_sel(i,6)<=30
        plot(EQ_sel(i,2),EQ_sel(i,5),'rx');
    elseif (EQ_sel(i,6)>30) && (EQ_sel(i,6)<=80)
        plot(EQ_sel(i,2),EQ_sel(i,5),'bx');
    else EQ_sel(i,6)>80
         plot(EQ_sel(i,2),EQ_sel(i,5),'gx');
    end
    
    if i==1
       hold on
    end
end
xlim([DOY_start DOY_end]);
ylim([0 max(EQ_sel(:,5))+10]);
%xlabel('Day of year');
ylabel('Seismic Index (K_{s})');
h=zeros(3,1);
h(1) = plot(NaN,NaN,'rx');
h(2) = plot(NaN,NaN,'bx');
h(3) = plot(NaN,NaN,'gx');
lgd=legend(h,'<30km','30-80km','>80km','location','northeast');
title(lgd,'Depth');
title(sprintf('Z/G Polarisation ratio | %s station, %d',station,year));
hold off

%Plotting polarisation ratio

subplot(5,1,2);
plot(days,ZG1,days,ZG1_mum2sig,'r--',days,ZG1_mup2sig,'r--');
%xlabel('Day of year');
ylabel('Z/G');
legend(sprintf('\\Deltaf=%.4f-%.4f Hz',f_A1,f_A2),'location','northeast');
xlim([DOY_start DOY_end]);

subplot(5,1,3);
plot(days,ZG2,days,ZG2_mum2sig,'r--',days,ZG2_mup2sig,'r--');
%xlabel('Day of year');
ylabel('Z/G');
legend(sprintf('\\Deltaf=%.4f-%.4f Hz',f_B1,f_B2),'location','northeast');
xlim([DOY_start DOY_end]);

subplot(5,1,4);
plot(days,ZG3,days,ZG3_mum2sig,'r--',days,ZG3_mup2sig,'r--');
%xlabel('Day of year');
ylabel('Z/G');
legend(sprintf('\\Deltaf=%.4f-%.4f Hz',f_C1,f_C2),'location','northeast');
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
