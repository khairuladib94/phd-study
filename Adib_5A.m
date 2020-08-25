stn1='CEB';                           %Abbreviation of station name
station1='Cebu';                       %Station full name
stn2='LGZ';
station2='Legazpi';

mag_min=5.0;                         %Minimum magnitude of earthquakes to be considered
dis_max=300;                         %Maximum epicentral distance from the station

%--------------------------------------------------------------------------

% datenum_start=round(min(UT1m));
% datenum_end=floor(max(UT1m));

date_start=[2011,12,8];           %Insert custom start and end dates
date_end=[2012,6,4];
datenum_start=datenum(date_start);
datenum_end=datenum(date_end);

%--------------------------------------------------------------------------\

days_num=datenum_end-datenum_start+1;
days_vec=datetime(datevec(datenum_start:1:datenum_end),'Format','dd/MM/yyyy');


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
stn_latlon1=[stn_IndPhi(stn_num1,2:3)];
stn_latlon2=[stn_IndPhi(stn_num2,2:3)];
stn_dis=deg2km(distance('gc',stn_latlon1,stn_latlon2));

%Earthquake----------------------------------------------------------------

%Building selected earthquakes table

j=1;
for i=1:length(EQ_IndPhi1117)
    if datenum(EQ_IndPhi1117(i,1:3))>=datenum_start && datenum(EQ_IndPhi1117(i,1:3))<=datenum_end
        EQ_table(j,:)=EQ_IndPhi1117(i,:);
        j=j+1;
    end
    if datenum(EQ_IndPhi1117(i,1:3))>datenum_end
        break
    end
end

i=1;
for j=1:length(EQ_table)
    if  EQ_table(j,10)>=mag_min
        EQ_latlon=[EQ_table(j,7:8)];
        EQ_dis=deg2km(distance('gc',stn_latlon1,EQ_latlon));
        EQ_date=datenum(EQ_table(j,1:3));
        EQ_mag=EQ_table(j,10);
        EQ_DOY=day(datetime([EQ_table(j,1:3)]),'dayofyear');
        EQ_Ks=((10^(0.75*EQ_mag))/(10*EQ_dis))*((1+EQ_dis*10^(-EQ_mag/2))^(-2.33));
        EQ_depth=EQ_table(j,9);
        EQ_f=1000/(pi*((EQ_depth*1000)^2)*1.2566e-06);
        
        if EQ_dis<=dis_max
            EQ_sel1(i,1)=EQ_date;
            EQ_sel1(i,2)=EQ_DOY;
            EQ_sel1(i,3)=EQ_dis;
            EQ_sel1(i,4)=EQ_mag;
            EQ_sel1(i,5)=EQ_Ks;
            EQ_sel1(i,6)=EQ_depth;
            EQ_sel1(i,7)=EQ_f;
            EQ_sel1(i,8:9)=EQ_latlon;
            i=i+1;
        end
    end
    
end


%%%%%%%%%%%%%


i=1;
for j=1:length(EQ_table)
    if  EQ_table(j,10)>=mag_min
        EQ_latlon=[EQ_table(j,7:8)];
        EQ_dis=deg2km(distance('gc',stn_latlon2,EQ_latlon));
        EQ_date=datenum(EQ_table(j,1:3));
        EQ_mag=EQ_table(j,10);
        EQ_DOY=day(datetime([EQ_table(j,1:3)]),'dayofyear');
        EQ_Ks=((10^(0.75*EQ_mag))/(10*EQ_dis))*((1+EQ_dis*10^(-EQ_mag/2))^(-2.33));
        EQ_depth=EQ_table(j,9);
        EQ_f=1000/(pi*((EQ_depth*1000)^2)*1.2566e-06);

        
        if EQ_dis<=dis_max
            EQ_sel2(i,1)=EQ_date;
            EQ_sel2(i,2)=EQ_DOY;
            EQ_sel2(i,3)=EQ_dis;
            EQ_sel2(i,4)=EQ_mag;
            EQ_sel2(i,5)=EQ_Ks;
            EQ_sel2(i,6)=EQ_depth;
            EQ_sel2(i,7)=EQ_f;
            EQ_sel2(i,8:9)=EQ_latlon;
            i=i+1;
        end
    end
    
end


%--------------------------------------------------------------------------

data_start=round(min(UT1m));
data_end=floor(max(UT1m));

day_start=datenum_start-data_start+1;
day_end=datenum_end-data_start+1;
datenum_vec=datenum_start:datenum_end;

H1_period=H1(86400*day_start-86399:86400*day_end);
D1_period=D1(86400*day_start-86399:86400*day_end);
Z1_period=Z1(86400*day_start-86399:86400*day_end);

H2_period=H2(86400*day_start-86399:86400*day_end);
D2_period=D2(86400*day_start-86399:86400*day_end);
Z2_period=Z2(86400*day_start-86399:86400*day_end);

%Removing noise/outlier
H1_dn=medfilt1(H1_period,'omitnan');
D1_dn=medfilt1(D1_period,'omitnan');
Z1_dn=medfilt1(Z1_period,'omitnan');

H2_dn=medfilt1(H2_period,'omitnan');
D2_dn=medfilt1(D2_period,'omitnan');
Z2_dn=medfilt1(Z2_period,'omitnan');

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

%Remove daily data which contain so many NaNs


%Calculating daily differences
for i=1:days_num
    
    H1_del(1,i)=max(H1_day(:,i))-min(H1_day(:,i));
    D1_del(1,i)=max(D1_day(:,i))-min(D1_day(:,i));
    Z1_del(1,i)=max(Z1_day(:,i))-min(Z1_day(:,i));
    
    H2_del(1,i)=max(H2_day(:,i))-min(H2_day(:,i));
    D2_del(1,i)=max(D2_day(:,i))-min(D2_day(:,i));
    Z2_del(1,i)=max(Z2_day(:,i))-min(Z2_day(:,i));

    j1=sum(isnan(H1_day(:,i)));
    k1=sum(isnan(D1_day(:,i)));
    l1=sum(isnan(Z1_day(:,i)));
    
    j2=sum(isnan(H2_day(:,i)));
    k2=sum(isnan(D2_day(:,i)));
    l2=sum(isnan(Z2_day(:,i)));
    
    if j1>50
        H1_del(1,i)=NaN;
    end
    if k1>50
        D1_del(1,i)=NaN;
    end
    if l1>50
        Z1_del(1,i)=NaN;
    end
    
    if j2>50
        H2_del(1,i)=NaN;
    end
    if k2>50
        D2_del(1,i)=NaN;
    end
    if l2>50
        Z2_del(1,i)=NaN;
    end

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
hold on
for i=1:size(EQ_sel1,1)
    j=datetime(datevec(EQ_sel1(i,1)),'Format','dd/MM/yyyy');
    plot(j,EQ_sel1(i,5),'ro');
end
for i=1:size(EQ_sel2,1)
    j=datetime(datevec(EQ_sel2(i,1)),'Format','dd/MM/yyyy');
    plot(j,EQ_sel2(i,5),'b^');
end
xlim([min(days_vec) max(days_vec)])
ylim([0 max(EQ_sel1(:,5))+10]);
ylabel('K_{s}');
set(gca,'Xticklabel',[]);
h=zeros(2,1);
h(1) = plot(days_vec,NaN(1,days_num),'ro');
h(2) = plot(days_vec,NaN(1,days_num),'b^');
lgd=legend(h,sprintf('%s',stn1),sprintf('%s',stn2),'location','northwest');
title(lgd,'Near to:');
title(sprintf('Diurnal variation | %s & %s station (%.0f km), %s - %s ',station1,station2,stn_dis,min(days_vec),max(days_vec)));
hold off

subplot(5,1,2);
plot(days_vec,H1_del,'r');
hold on
plot(days_vec,H2_del,'b--');
hold off
xlim([min(days_vec) max(days_vec)])
legend('CEB','LGZ','location','northwest');
set(gca,'Xticklabel',[]);
ylabel('\DeltaH (nT)');
title(sprintf('Diurnal variation | %s & %s station (%.0f km), %s - %s ',station1,station2,stn_dis,min(days_vec),max(days_vec)));


subplot(5,1,3);
plot(days_vec,D1_del,'r');
hold on
plot(days_vec,D2_del,'b--');
hold off
xlim([min(days_vec) max(days_vec)])
legend('CEB','LGZ','location','northwest');
set(gca,'Xticklabel',[]);
ylabel('\DeltaD (nT)');


subplot(5,1,4);
plot(days_vec,Z1_del,'r');
hold on
plot(days_vec,Z2_del,'b--');
hold off
xlim([min(days_vec) max(days_vec)])
legend('CEB','LGZ','location','northwest');
ylabel('\DeltaZ (nT)');
xlabel('Date');

%Getting Kp and Dst indices
j=1;
for i=1:length(index_1117)
    m=datenum(index_1117(i,1),1,index_1117(i,2));
    if m>=datenum_start && m<=datenum_end
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
bar(days_vec,Kp,'b');
hold on
bar(days_vec,Dst,'r');
plot(days_vec,ones(days_num,1)*(32),'b--');
plot(days_vec,ones(days_num,1)*(-50),'r--');
hold off
xlabel('Date');
ylabel('Dst     \SigmaK_{p}');
legend('\SigmaK_{p}','Dst','location','southwest');
xlim([min(days_vec) max(days_vec)])
ylim([-60 40]);

%------------



% subplot(1,2,1)
% plot(min_vec,H1_day(:,40));
% datetick('x','HH:MM');
% title('CEB | 16/1/2012');
% ylabel('H (nT)');
% xlabel('Local time');
% 
% subplot(1,2,2)
% plot(min_vec,H1_day(:,41));
% datetick('x','HH:MM');
% title('CEB | 17/1/2012');
% ylabel('H (nT)');
% xlabel('Local time');

