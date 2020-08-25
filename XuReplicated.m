stn1='CEB';                           %Abbreviation of station name
station1='Cebu';                %Station full name
stn2='LGZ';
station2='Legazpi';

mag_min=5.0;                         %Minimum magnitude of earthquakes to be considered
dis_max=300;                         %Maximum epicentral distance from the station

%%

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

%Setting station number
if stn1=='TGG'
    stn1_num=1;
elseif stn1=='MUT'
    stn_num1=2;
elseif stn1=='LGZ'
    stn_num1=3;
elseif stn1=='CEB'
    stn_num1=4;
elseif stn1=='CDO'
    stn_num1=5;
elseif stn1=='DAV'
    stn_num1=6;
elseif stn1=='GSI'
    stn_num1=7;
elseif stn1=='SCN'
    stn_num1=8;
elseif stn1=='LWA'
    stn_num1=9;
elseif stn1=='PTN'
    stn_num1=10;
elseif stn1=='MND'
    stn_num1=11;
elseif stn1=='BIK'
    stn_num1=12;
elseif stn1=='JYP'
    stn_num1=13;
elseif stn1=='PRP'
    stn_num1=14;
elseif stn1=='KPG'
    stn_num1=15;
elseif stn1=='KTB'
    stn_num1=16;
elseif stn1=='DAW'
    stn_num1=17;
end

if stn2=='TGG'
    stn2_num=1;
elseif stn2=='MUT'
    stn_num2=2;
elseif stn2=='LGZ'
    stn_num2=3;
elseif stn2=='CEB'
    stn_num2=4;
elseif stn2=='CDO'
    stn_num2=5;
elseif stn2=='DAV'
    stn_num2=6;
elseif stn2=='GSI'
    stn_num2=7;
elseif stn2=='SCN'
    stn_num2=8;
elseif stn2=='LWA'
    stn_num2=9;
elseif stn2=='PTN'
    stn_num2=10;
elseif stn2=='MND'
    stn_num2=11;
elseif stn2=='BIK'
    stn_num2=12;
elseif stn2=='JYP'
    stn_num2=13;
elseif stn2=='PRP'
    stn_num2=14;
elseif stn2=='KPG'
    stn_num2=15;
elseif stn2=='KTB'
    stn_num2=16;
elseif stn2=='DAW'
    stn_num2=17;
end

%Getting station coordinate
stn_latlon1=[stn_IndPhi(stn_num1,2) stn_IndPhi(stn_num1,3)];
stn_latlon2=[stn_IndPhi(stn_num2,2) stn_IndPhi(stn_num2,3)];
stn_dis=deg2km(distance('gc',stn_latlon1,stn_latlon2));

%%

%Building selected earthquakes table
EQ_sel1=NaN(100,7);
i=1;
for j=1:length(EQ_IndPhi1117)
    if EQ_IndPhi1117(j,1)==year
        if  EQ_IndPhi1117(j,10)>=mag_min
            EQ_latlon=[EQ_IndPhi1117(j,7) EQ_IndPhi1117(j,8)];
            EQ_dis=deg2km(distance('gc',stn_latlon1,EQ_latlon));
            EQ_date=datenum([EQ_IndPhi1117(j,1) EQ_IndPhi1117(j,2) EQ_IndPhi1117(j,3)]);
            EQ_mag=EQ_IndPhi1117(j,10);
            EQ_DOY=day(datetime([EQ_IndPhi1117(j,1) EQ_IndPhi1117(j,2) EQ_IndPhi1117(j,3)]),'dayofyear');
            EQ_Ks=((10^(0.75*EQ_mag))/(10*EQ_dis))*((1+EQ_dis*10^(-EQ_mag/2))^(-2.33));
            EQ_depth=EQ_IndPhi1117(j,9);
            EQ_f=(2*600)/(((1000*EQ_depth)^2)*1.2566e-6*2*pi);
            
            if EQ_dis<=dis_max
                EQ_sel1(i,1)=EQ_date;
                EQ_sel1(i,2)=EQ_DOY;
                EQ_sel1(i,3)=EQ_dis;
                EQ_sel1(i,4)=EQ_mag;
                EQ_sel1(i,5)=EQ_Ks;
                EQ_sel1(i,6)=EQ_depth;
                EQ_sel1(i,7)=EQ_f;
                i=i+1;
            end
        end
    end
end

EQ_sel1(any(isnan(EQ_sel1),2),:)=[];
if size(EQ_sel1,1)>1
    EQ_sel1=flip(EQ_sel1);
end

%%%

EQ_sel2=NaN(100,7);
i=1;
for j=1:length(EQ_IndPhi1117)
    if EQ_IndPhi1117(j,1)==year
        if  EQ_IndPhi1117(j,10)>=mag_min
            EQ_latlon=[EQ_IndPhi1117(j,7) EQ_IndPhi1117(j,8)];
            EQ_dis=deg2km(distance('gc',stn_latlon2,EQ_latlon));
            EQ_date=datenum([EQ_IndPhi1117(j,1) EQ_IndPhi1117(j,2) EQ_IndPhi1117(j,3)]);
            EQ_mag=EQ_IndPhi1117(j,10);
            EQ_DOY=day(datetime([EQ_IndPhi1117(j,1) EQ_IndPhi1117(j,2) EQ_IndPhi1117(j,3)]),'dayofyear');
            EQ_Ks=((10^(0.75*EQ_mag))/(10*EQ_dis))*((1+EQ_dis*10^(-EQ_mag/2))^(-2.33));
            EQ_depth=EQ_IndPhi1117(j,9);
            EQ_f=(2*600)/(((1000*EQ_depth)^2)*1.2566e-6*2*pi);
            
            if EQ_dis<=dis_max
                EQ_sel2(i,1)=EQ_date;
                EQ_sel2(i,2)=EQ_DOY;
                EQ_sel2(i,3)=EQ_dis;
                EQ_sel2(i,4)=EQ_mag;
                EQ_sel2(i,5)=EQ_Ks;
                EQ_sel2(i,6)=EQ_depth;
                EQ_sel2(i,7)=EQ_f;
                i=i+1;
            end
        end
    end
end

EQ_sel2(any(isnan(EQ_sel2),2),:)=[];
if size(EQ_sel2,1)>1
    EQ_sel2=flip(EQ_sel2);
end

%%

%Removing noise/outlier
H1_dn=medfilt1(H1,'omitnan');
D1_dn=medfilt1(D1,'omitnan');
Z1_dn=medfilt1(Z1,'omitnan');

H2_dn=medfilt1(H2,'omitnan');
D2_dn=medfilt1(D2,'omitnan');
Z2_dn=medfilt1(Z2,'omitnan');

for i=1:5
    
    H1_sig=std(H1_dn,'omitnan');
    H1_mu=mean(H1_dn,'omitnan');
    D1_sig=std(D1_dn,'omitnan');
    D1_mu=mean(D1_dn,'omitnan');
    Z1_sig=std(Z1_dn,'omitnan');
    Z1_mu=mean(Z1_dn,'omitnan');
    
    H2_sig=std(H2_dn,'omitnan');
    H2_mu=mean(H2_dn,'omitnan');
    D2_sig=std(D2_dn,'omitnan');
    D2_mu=mean(D2_dn,'omitnan');
    Z2_sig=std(Z2_dn,'omitnan');
    Z2_mu=mean(Z2_dn,'omitnan');
    
    for j=1:length(H1_dn)
        
        if H1_dn(j)>H1_mu+5*H1_sig||H1_dn(j)<H1_mu-5*H1_sig
            H1_dn(j)=NaN;
        end
        if D1_dn(j)>D1_mu+5*D1_sig||D1_dn(j)<D1_mu-5*D1_sig
            D1_dn(j)=NaN;
        end
        if Z1_dn(j)>Z1_mu+5*Z1_sig||Z1_dn(j)<Z1_mu-5*Z1_sig
            Z1_dn(j)=NaN;
        end
        
        if H2_dn(j)>H2_mu+5*H2_sig||H2_dn(j)<H2_mu-5*H2_sig
            H2_dn(j)=NaN;
        end
        if D2_dn(j)>D2_mu+5*D2_sig||D2_dn(j)<D2_mu-5*D2_sig
            D2_dn(j)=NaN;
        end
        if Z2_dn(j)>Z2_mu+5*Z2_sig||Z2_dn(j)<Z2_mu-5*Z2_sig
            Z2_dn(j)=NaN;
        end
        
    end
end

%%

%Downsample to minute data
H1_ds=downsample(H1_dn,60);
D1_ds=downsample(D1_dn,60);
Z1_ds=downsample(Z1_dn,60);

H2_ds=downsample(H2_dn,60);
D2_ds=downsample(D2_dn,60);
Z2_ds=downsample(Z2_dn,60);

%Local time zone correction
t_local=7;
t_beyond=24-t_local;
t_start=t_beyond*60+1;
t_end=length(H1_ds)-t_local*60;

H1_lt=NaN(1,days_num*1440);
D1_lt=NaN(1,days_num*1440);
Z1_lt=NaN(1,days_num*1440);

% H2_lt=NaN(days_num*1440,1);
D2_lt=NaN(days_num*1440,1);
Z2_lt=NaN(days_num*1440,1);

H1_lt(1,1441:days_num*1440)=H1_ds(t_start:t_end);
D1_lt(1,1441:days_num*1440)=D1_ds(t_start:t_end);
Z1_lt(1,1441:days_num*1440)=Z1_ds(t_start:t_end);

H2_lt(1441:days_num*1440)=H2_ds(t_start:t_end);
D2_lt(1441:days_num*1440)=D2_ds(t_start:t_end);
Z2_lt(1441:days_num*1440)=Z2_ds(t_start:t_end);

%Reshaping into daily data
H1_day=reshape(H1_lt,1440,[]);
D1_day=reshape(D1_lt,1440,[]);
Z1_day=reshape(Z1_lt,1440,[]);

H2_day=reshape(H2_lt,1440,[]);
D2_day=reshape(D2_lt,1440,[]);
Z2_day=reshape(Z2_lt,1440,[]);

%Calculating daily differences
H1_del=NaN(1,days_num);
D1_del=NaN(1,days_num);
Z1_del=NaN(1,days_num);

H2_del=NaN(1,days_num);
D2_del=NaN(1,days_num);
Z2_del=NaN(1,days_num);

for i=1:days_num
    
    H1_del(1,i)=max(H1_day(:,i))-min(H1_day(:,i));
    D1_del(1,i)=max(D1_day(:,i))-min(D1_day(:,i));
    Z1_del(1,i)=max(Z1_day(:,i))-min(Z1_day(:,i));
    
    H2_del(1,i)=max(H2_day(:,i))-min(H2_day(:,i));
    D2_del(1,i)=max(D2_day(:,i))-min(D2_day(:,i));
    Z2_del(1,i)=max(Z2_day(:,i))-min(Z2_day(:,i));
    
end


%Ratio of station 1:station 2
H1H2=H1_del./H2_del;
D1D2=D1_del./D2_del;
Z1Z2=Z1_del./Z2_del;

%mu+/-3sigma
H1H2_mu=mean(H1H2,'omitnan');
D1D2_mu=mean(D1D2,'omitnan');
Z1Z2_mu=mean(Z1Z2,'omitnan');

H1H2_sig=std(H1H2,'omitnan');
D1D2_sig=std(D1D2,'omitnan');
Z1Z2_sig=std(Z1Z2,'omitnan');

H1H2_mup3sig=H1H2_mu+3*H1H2_sig;
D1D2_mup3sig=D1D2_mu+3*D1D2_sig;
Z1Z2_mup3sig=Z1Z2_mu+3*Z1Z2_sig;

H1H2_mum3sig=H1H2_mu-3*H1H2_sig;
D1D2_mum3sig=D1D2_mu-3*D1D2_sig;
Z1Z2_mum3sig=Z1Z2_mu-3*Z1Z2_sig;

%----------------------------------------------------------------------------------------------------------------------------------------------
%Plotting


%Plotting seismic index chart
subplot(5,1,1)
for i=1:size(EQ_sel1,1)
    if EQ_sel1(i,6)<=30
        plot(EQ_sel1(i,2),EQ_sel1(i,5),'rx');
    elseif (EQ_sel1(i,6)>30) && (EQ_sel1(i,6)<=80)
        plot(EQ_sel1(i,2),EQ_sel1(i,5),'bx');
    else EQ_sel1(i,6)>80
        plot(EQ_sel1(i,2),EQ_sel1(i,5),'gx');
    end
    
    if i==1
        hold on
    end
end
for i=1:size(EQ_sel2,1)
    if EQ_sel2(i,6)<=30
        plot(EQ_sel2(i,2),EQ_sel2(i,5),'ro');
    elseif (EQ_sel2(i,6)>30) && (EQ_sel2(i,6)<=80)
        plot(EQ_sel2(i,2),EQ_sel2(i,5),'bo');
    else EQ_sel2(i,6)>80
        plot(EQ_sel2(i,2),EQ_sel2(i,5),'go');
    end
end
xlim([DOY_start DOY_end]);
ylim([0 max(EQ_sel1(:,5))+10]);
ylabel('Seismic Index (K_{s})');
h=zeros(2,1);
h(1) = plot(NaN,NaN,'kx');
h(2) = plot(NaN,NaN,'ko');
lgd=legend(h,sprintf('%s',stn1),sprintf('%s',stn2),'location','northeast');
title(lgd,'Near to:');
title(sprintf('Diurnal variation ratio | %d | %s-%s station (%.0f km)',year,station1,station2,stn_dis));
hold off



subplot(5,1,2);
plot(days,H1H2,days,ones(1,days_num)*H1H2_mup3sig,'r--',days,ones(1,days_num)*H1H2_mum3sig,'r--');
xlim([1 days_num]);
ylabel(sprintf('\\DeltaH_{%s}/\\DeltaH_{%s}',stn1,stn2));

subplot(5,1,3);
plot(days,D1D2,days,ones(1,days_num)*D1D2_mup3sig,'r--',days,ones(1,days_num)*D1D2_mum3sig,'r--');
xlim([1 days_num]);
ylabel(sprintf('\\DeltaD_{%s}/\\DeltaD_{%s}',stn1,stn2));


subplot(5,1,4);
plot(days,Z1Z2,days,ones(1,days_num)*Z1Z2_mup3sig,'r--',days,ones(1,days_num)*Z1Z2_mum3sig,'r--');
xlim([1 days_num]);
ylabel(sprintf('\\DeltaZ_{%s}/\\DeltaZ_{%s}',stn1,stn2));

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



