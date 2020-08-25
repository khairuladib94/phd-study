mag_min=5.5;
dis_max=500;

date_start=[2007,1,1];           
date_end=[2018,3,23];
datenum_start=datenum(date_start);
datenum_end=datenum(date_end);

days_num=datenum_end-datenum_start+1;
days_vec=datetime(datevec(datenum_start:1:datenum_end),'Format','dd/MM/yyyy');

j=1;
for i=1:length(EQ_WORLD)
    if datenum(EQ_WORLD(i,1:3))>=datenum_start && datenum(EQ_WORLD(i,1:3))<=datenum_end
        EQ_table(j,:)=EQ_WORLD(i,:);
        j=j+1;
    end
    if datenum(EQ_WORLD(i,1:3))>datenum_end
        break
    end
end
    
    

i=1;
for j=1:length(EQ_table)
    if  EQ_table(j,10)>=mag_min
        EQ_latlon=[EQ_table(j,7:8)];
        EQ_dis=10000;
        for k=1:size(stn_MAGDAS,1)
            stn_latlon=[stn_MAGDAS(k,2:3)];
            EQ_dis_test=deg2km(distance('gc',stn_latlon,EQ_latlon));
            if EQ_dis_test<=EQ_dis
                EQ_dis=EQ_dis_test;
            end
        end
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

worldmap([-90 90],[-180 180])
geoshow('landareas.shp','FaceColor','White')
title(sprintf('World Earthquake Map | %s - %s, M > %0.1f, d < %.0f km',min(days_vec),max(days_vec),mag_min,dis_max));

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
    date_text=datetime(datevec(EQ_sel(i,1)),'Format','dd/MM/yyyy');
    if (EQ_sel(i,4)>=6.0 && EQ_sel(i,3)<=100) || EQ_sel(i,3)<=30 || EQ_sel(i,4)>=7.0
        textm(EQ_lat,EQ_lon,sprintf('%s(%0.1f)',date_text,EQ_sel(i,4)),'FontSize',8);
    end
end

for s=1:size(stn_MAGDAS,1)
    
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
    elseif s==22
        stn='KTN';
    elseif s==23
        stn='TIK';
    elseif s==24
        stn='CHD';
    elseif s==25
        stn='CST';
    elseif s==26
        stn='ZYK';
    elseif s==27
        stn='ZGN';
    elseif s==28
        stn='MGD';
    elseif s==29
        stn='YAK';
    elseif s==30
        stn='PTK';
    elseif s==31
        stn='ASB';
    elseif s==32
        stn='TNO';
    elseif s==33
        stn='ONW';
    elseif s==34
        stn='OIS';
    elseif s==35
        stn='KUJ';
    elseif s==36
        stn='AMA';
    elseif s==37
        stn='HLN';
    elseif s==38
        stn='EWA';
    elseif s==39
        stn='YAP';
    elseif s==40
        stn='BCL';
    elseif s==41
        stn='HVD';
    elseif s==42
        stn='TIR';
    elseif s==43
        stn='CMB';
    elseif s==44
        stn='CKT';
    elseif s==45
        stn='TWV';
    elseif s==46
        stn='ROC';
    elseif s==47
        stn='LMT';
    elseif s==48
        stn='CGR';
    elseif s==49
        stn='CMD';
    elseif s==50
        stn='CAN';
    elseif s==51
        stn='MLB';
    elseif s==52
        stn='HOB';
    elseif s==53
        stn='MCQ';
    elseif s==54
        stn='DVS';
    elseif s==55
        stn='WAD';
    elseif s==56
        stn='GLY';
    elseif s==57
        stn='JRS';
    elseif s==58
        stn='TPT';
    elseif s==59
        stn='TMA';
    elseif s==60
        stn='ANC';
    elseif s==61
        stn='HUA';
    elseif s==62
        stn='ICA';
    elseif s==63
        stn='EUS';
    elseif s==64
        stn='SMA';
    elseif s==65
        stn='LAQ';
    elseif s==66
        stn='FYM';
    elseif s==67
        stn='ASW';
    elseif s==68
        stn='KRT';
    elseif s==69
        stn='AAB';
    elseif s==70
        stn='ILR';
    elseif s==71
        stn='ABU';
    elseif s==72
        stn='LAG';
    elseif s==73
        stn='ABJ';
    elseif s==74
        stn='NAB';
    elseif s==75
        stn='DES';
    elseif s==76
        stn='LSK';
    elseif s==77
        stn='MPT';
    elseif s==78
        stn='DRB';
    elseif s==79
        stn='HER';
    end
    
    stn_latlon=[stn_MAGDAS(s,2:3)];
    plotm(stn_latlon(1),stn_latlon(2),'k^','MarkerFaceColor','k','MarkerSize',5);
    textm(stn_latlon(1),stn_latlon(2)+0.2,stn)
 
end

