DOY_start=195;
DOY_end=205;

year=min(datevec(UT1m));
year=year(1,1);

days_num=DOY_end-DOY_start+1;
days=DOY_start:1:DOY_end;

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


H_welch=NaN(451,days_num);
i=1;
f_s=1;
for i=1:days_num
    [H_welch(:,i),f]=pwelch(H_night(:,i),1800,[],900,f_s);
end
    
plot(f,10*log10(H_welch));

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

ZH1=Z_norm./H_norm;

%mu+/-2sigma
ZH1_mu=mean(ZH1,'omitnan');
ZH1_sig=std(ZH1,'omitnan');
ZH1_mup2sig=ones(1,days_num)*(ZH1_mu+2*ZH1_sig);
ZH1_mum2sig=ones(1,days_num)*(ZH1_mu-2*ZH1_sig);
ZH1_mm=movmean(ZH1,5,'omitnan');