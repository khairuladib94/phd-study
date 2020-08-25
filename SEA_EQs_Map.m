mag_min=7.0;
dis_max=600;

date_start=[2010,1,1];           %Insert custom start and end dates
date_end=[2018,3,14];
datenum_start=datenum(date_start);
datenum_end=datenum(date_end);

days_num=datenum_end-datenum_start+1;
days_vec=datetime(datevec(datenum_start:1:datenum_end),'Format','dd/MM/yyyy');

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

for k=1:17
    
    stn_latlon=[stn_IndPhi(k,2:3)];

    i=1;
    for j=1:length(EQ_table)
        if  EQ_table(j,10)>=mag_min
            EQ_latlon=[EQ_table(j,7:8)];
            EQ_dis=deg2km(distance('gc',stn_latlon,EQ_latlon));
            EQ_date=datenum(EQ_table(j,1:3));
            EQ_mag=EQ_table(j,10);
            EQ_DOY=day(datetime([EQ_table(j,1:3)]),'dayofyear');
            EQ_Ks=((10^(0.75*EQ_mag))/(10*EQ_dis))*((1+EQ_dis*10^(-EQ_mag/2))^(-2.33));
            EQ_depth=EQ_table(j,9);
            EQ_f=1000/(pi*((EQ_depth*1000)^2)*1.2566e-06);
            
            if EQ_dis<=dis_max
                EQ_sel(i,1)=EQ_date;
                EQ_sel(i,2)=EQ_DOY;
                EQ_sel(i,3)=EQ_dis;
                EQ_sel(i,4)=EQ_mag;
                EQ_sel(i,5)=EQ_Ks;
                EQ_sel(i,6)=EQ_depth;
                EQ_sel(i,7)=EQ_f;
                EQ_sel(i,8:9)=EQ_latlon;
                i=i+1;
            end
            
        end
        
    end
    
end

worldmap([-10 20],[95 128])
geoshow('landareas.shp','FaceColor','White')

for i=1:size(EQ_sel,1)
    EQ_lat=EQ_sel(i,8);
    EQ_lon=EQ_sel(i,9);
    EQ_radius=((EQ_sel(i,4)-4.5)/0.025);
    if EQ_sel(i,6)<=50
        colour='r';
    elseif (EQ_sel(i,6)>50) && (EQ_sel(i,6)<=200)
        colour='b';
    else EQ_sel(i,6)>200
         colour='g';
    end
    circlem(EQ_lat,EQ_lon,EQ_radius,'edgecolor','none','facecolor',colour,'facealpha','0.2');
end

for s=1:17
    
    if  s==1
        stn='TGG';
    elseif s==2
        stn='MUT';
    elseif s==3
        stn='LGZ';
    elseif s==4
        stn='CEB';
    elseif s==5
        stn='CDO';
    elseif s==6
        stn='DAV';
    elseif s==7
        stn='GSI';
    elseif s==8
        stn='SCN';
    elseif s==9
        stn='LWA';
    elseif s==10
        stn='PTN';
    elseif s==11 
        stn='MND';
    elseif s==12
        stn='BIK';
    elseif s==13
        stn='JYP';
    elseif s==14
        stn='PRP';
    elseif s==15
        stn='KPG';
    elseif s==16
        stn='KTB';
    elseif s==17
        stn='DAW';
    elseif s==18
        stn='LKW';
    elseif s==19
        stn='SBH';
    elseif s==20
        stn='PER';
    elseif s==21
        stn='BTG';
    end
   
    stn_latlon=[stn_IndPhi(s,2:3)];
    plotm(stn_latlon(1),stn_latlon(2),'k^','MarkerFaceColor','k','MarkerSize',5);
 
end

