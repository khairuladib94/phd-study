year=2017;
mag_min=5.0;                         
dis_max=300;

DOY_start=1;
if rem(year,4)==0
    DOY_end=366;
else
    DOY_end=365;
end
days_num=DOY_end;
days=DOY_start:1:DOY_end;

for i=1:16
    
    if i==1
        stn='TGG';
    elseif i==2
        stn='MUT';
    elseif i==3
        stn='LGZ';
    elseif i==4
        stn='CEB';
    elseif i==5
        stn='CDO';
    elseif i==6
        stn='DAV';
    elseif i==7;
        stn='GSI';
    elseif i==8
        stn='SCN';
    elseif i==9
        stn='LWA';
    elseif i==10
        stn='PTN';
    elseif i==11
        stn='MND';
    elseif i==12
        stn='BIK';
    elseif i==13
        stn='JYP';
    elseif i==14
        stn='PRP';
    elseif i==15
        stn='KPG';
    elseif i==16
        stn='KTB';
    end
    
    stn_latlon=[stn_IndPhi(i,2) stn_IndPhi(i,3)];
    
    EQ_sel=NaN(100,7);
    k=1;
    for l=1:length(EQ_IndPhi1117)
        if EQ_IndPhi1117(l,1)==year
            if  EQ_IndPhi1117(l,10)>=mag_min
                EQ_latlon=[EQ_IndPhi1117(l,7) EQ_IndPhi1117(l,8)];
                EQ_dis=deg2km(distance('gc',stn_latlon,EQ_latlon));
                EQ_date=datenum([EQ_IndPhi1117(l,1) EQ_IndPhi1117(l,2) EQ_IndPhi1117(l,3)]);
                EQ_mag=EQ_IndPhi1117(l,10);
                EQ_DOY=day(datetime([EQ_IndPhi1117(l,1) EQ_IndPhi1117(l,2) EQ_IndPhi1117(l,3)]),'dayofyear');
                EQ_Ks=((10^(0.75*EQ_mag))/(10*EQ_dis))*((1+EQ_dis*10^(-EQ_mag/2))^(-2.33));
                EQ_depth=EQ_IndPhi1117(l,9);
                EQ_f=(2*600)/(((1000*EQ_depth)^2)*1.2566e-6*2*pi);
                
                if EQ_dis<=dis_max
                    EQ_sel(k,1)=EQ_date;
                    EQ_sel(k,2)=EQ_DOY;
                    EQ_sel(k,3)=EQ_dis;
                    EQ_sel(k,4)=EQ_mag;
                    EQ_sel(k,5)=EQ_Ks;
                    EQ_sel(k,6)=EQ_depth;
                    EQ_sel(k,7)=EQ_f;
                    k=k+1;
                end
            end
        end
    end
    
    EQ_sel(any(isnan(EQ_sel),2),:)=[];
    if size(EQ_sel,1)>1
        EQ_sel=flip(EQ_sel);
    end
    
    if i<=8
        if i==1
            figure(1)
        end
        subplot(4,2,i)
    end
    
    if i>8
        if i==9
            figure(2)
        end
        j=i-8;
        subplot(4,2,j)
    end
    
    for m=1:size(EQ_sel,1)
        if EQ_sel(m,6)<=30
            plot(EQ_sel(m,2),EQ_sel(m,5),'rx');
        elseif (EQ_sel(m,6)>30) && (EQ_sel(m,6)<=80)
            plot(EQ_sel(m,2),EQ_sel(m,5),'bx');
        else EQ_sel(m,6)>80
            plot(EQ_sel(m,2),EQ_sel(m,5),'gx');
        end
        if m==1
            hold on
        end
    end
    xlim([DOY_start DOY_end]);
    title(sprintf('%s',stn));
    hold off
        
end
    
figure(1)
suptitle(sprintf('%d',year));
figure(2)
suptitle(sprintf('%d',year));


        
    