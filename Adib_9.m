%Place and time inputs-----------------------------------------------------

stn='TNO';                       %Abbreviation of station name
date_start=[2015,1,1];           %Insert custom start and end dates
date_end=[2015,12,31];           %Period spanning through 3 consecutive years is the maximum  

%Customization-------------------------------------------------------------

mag_min=5.0;                     %Minimum magnitude of earthquakes to be considered
dis_max=300;                     %Maximum epicentral distance from the station

f_A1=0.01;
f_A2=0.1;

LT_start=0;
LT_end=4;

%Time period---------------------------------------------------------------

datenum_start=datenum(date_start);
datenum_end=datenum(date_end);
if date_start(1)==date_end(1)
    year=string(date_start(1));
else
    year=sprintf('%s-%s',string(date_start(1)),string(date_end(1)));
end

%Loading files-------------------------------------------------------------

load VARIABLES_WORLD

j=date_start(1);
for i=1:3
    year_vec(1,i)=j;
    if j==date_end(1)
        break
    end
    j=j+1;
end

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

%--------------------------------------------------------------------------

days_num=datenum_end-datenum_start+1;
days_vec=datetime(datevec(datenum_start:1:datenum_end),'Format','dd/MM/yyyy');

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


%Geomagnetic data----------------------------------------------------------

data_start=round(min(UT1m));
data_end=floor(max(UT1m));

day_start=datenum_start-data_start+1;
day_end=datenum_end-data_start+1;
datenum_vec=datenum_start:datenum_end;

H_period=H(86400*day_start-86399:86400*day_end);
D_period=D(86400*day_start-86399:86400*day_end);
Z_period=Z(86400*day_start-86399:86400*day_end);

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

%Nighttime data in LT
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
t_starts=t_start*3600;      
t_int=LT_length*3600;

for i=1:days_num
    
    H_night(:,i)=H_dn(t_starts:t_starts+t_int-1);
    D_night(:,i)=D_dn(t_starts:t_starts+t_int-1);
    Z_night(:,i)=Z_dn(t_starts:t_starts+t_int-1);
    t_starts=t_starts+86400;
    
end

%---------------------------------------------------------

ano1=300;
ano2=350;

% for i=ano1:ano2
%     
%     H_filt=filter1('bp',H_night(:,i),'fc',[0.01 0.1],'fs',1);
%     D_filt=filter1('bp',D_night(:,i),'fc',[0.01 0.1],'fs',1);
%     
%     [H_up,H_down]=envelope(H_filt,100,'peak');
%     [D_up,D_down]=envelope(D_filt,100,'peak');
%     
%     theta(:,i)=rad2deg(atan(H_up./D_up));
%     
% end
 


H_hil=hilbert(H_night);
D_hil=hilbert(D_night);

for i=1:size(H_hil,2)
    
    theta(:,i)=(atan(((2*abs(H_hil(:,i)).*abs(D_hil(:,i)))./((abs(D_hil(:,i)).^2)-(abs(H_hil(:,i)).^2))).*cos(angle(H_hil(:,i))-angle(D_hil(:,i)))))/2;
    
end

plot(theta)




for i=1:length(theta)
    if H_nightmu(1,i)>0
        if D_nightmu(1,i)>0
            alpha(1,i)=theta(1,i);
        else
            alpha(1,i)=pi-theta(1,i);
        end
    else     
        if D_nightmu(1,i)>0
            alpha(1,i)=2*pi-theta(1,i);
        else
            alpha(1,i)=pi+theta(1,i);
        end 
    end
end



%---------------------------------



%Plotting phase angles
subplot(4,1,4)
plot(days_vec,alpha)
xlabel('Date')
ylabel('Azimuthal angle (\theta)')

set(gcf, 'Position', get(0, 'Screensize'));

%--------------------------------------------------------------------------
