year=2013;
stn='DAV';

if stn=='TGG'
    stnnum=1;
end
if stn=='MUT'
    stnnum=2;
end
if stn=='LGZ'
    stnnum=3;
end
if stn=='CEB'
    stnnum=4;
end
if stn=='CDO'
    stnnum=5;
end
if stn=='DAV'
    stnnum=6;
end
if stn=='GSI'
    stnnum=7;
end
if stn=='SCN'
    stnnum=8;
end
if stn=='LWA'
    stnnum=9;
end
if stn=='PTN'
    stnnum=10;
end
if stn=='MND';
stnnum=11;
end
if stn=='BIK'
    stnnum=12;
end
if stn=='JYP'
    stnnum=13;
end
if stn=='PRP'
    stnnum=14;
end
if stn=='KPG'
    stnnum=15;
end
    
latlonstn=[StnIndPhillip(stnnum,2) StnIndPhillip(stnnum,3)];
    
%check leap year
if rem(year,4)==0
    days=1:1:366;
    numdays=366;
else
    days=1:1:365;
    numdays=365;
end

Eq20112017=vertcat(EqIndo20112017,EqPhillip20112017);
Eq20112017=sortrows(Eq20112017,3);
Eq20112017=sortrows(Eq20112017,2);
Eq20112017=sortrows(Eq20112017,1);

EqInd=zeros(numdays,1);
EqInd2=zeros(numdays,1);
k=1;

for i=1:length(Eq20112017)
    
    if (Eq20112017(i,1)==year)
        
        daynum=day(datetime([Eq20112017(i,1) Eq20112017(i,2) Eq20112017(i,3)]),'dayofyear');
        latloneq=[Eq20112017(i,7) Eq20112017(i,8)];
        dis=deg2km(distance('gc',latlonstn,latloneq));
        mag=Eq20112017(i,10);
        idx=mag;
        
        if (dis<=5000) && (daynum==k)
            EqInd(k)=EqInd(k)+idx;
            
        end
        if (dis<=5000) && (daynum~=k)
            k=daynum;
            EqInd(k)=EqInd(k)+idx;
        end
    
    end
    
end
    
    
%plot(days,EqInd,days,EqInd2,'black');

hold on

for i=1:length(Eq20112017)
    
    if (Eq20112017(i,1)==year)
        
        daynum=day(datetime([Eq20112017(i,1) Eq20112017(i,2) Eq20112017(i,3)]),'dayofyear');
        latloneq=[Eq20112017(i,7) Eq20112017(i,8)];
        dis=deg2km(distance('gc',latlonstn,latloneq));
        mag=Eq20112017(i,10);
        
        if dis<=5000
            
            if mag<3.0
                plot(daynum,mag,'green*');
                hold on
            end
            
            if mag>=3.0 && mag<6.0
                plot(daynum,mag,'yellow*');
                hold on
            end
            
            if mag>=6.0
                plot(daynum,mag,'red*');
                hold on
            end
            
            
        end
        
    
    end
    
end
hold off

axis([1 numdays -inf inf]);
            
            
            
            
     
