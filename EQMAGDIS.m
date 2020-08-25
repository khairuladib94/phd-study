Eq20112017=vertcat(EqIndo20112017,EqPhillip20112017);
ctsize=length(Eq20112017)*length(StnIndPhillip);
comptable=string(zeros(ctsize,5));

i=1;
for j=1:length(StnIndPhillip)
     
    for k=1:length(Eq20112017)
        
        if j==1
            stn='TGG';
        end
        if j==2
            stn='MUT';
        end
        if j==3
            stn='LGZ';
        end
        if j==4
            stn='CEB';
        end
        if j==5
            stn='CDO';
        end
        if j==6
            stn='DAV';
        end
        if j==7
            stn='GSI';
        end
        if j==8
            stn='SCN';
        end
        if j==9
            stn='LWA';
        end
        if j==10
            stn='PTN';
        end
        if j==11
            stn='MND';
        end
        if j==12
            stn='BIK';
        end
        if j==13
            stn='JYP';
        end
        if j==14
            stn='PRP';
        end
        if j==15
            stn='KPG';
        end
        
        dateno=datenum([Eq20112017(k,1) Eq20112017(k,2) Eq20112017(k,3)]);
        date=datetime(datevec(dateno),'Format','dd/MM/yyy');
        latlonstn=[StnIndPhillip(j,2) StnIndPhillip(j,3)];
        latloneq=[Eq20112017(k,7) Eq20112017(k,8)];
        dis=deg2km(distance('gc',latlonstn,latloneq));
        mag=Eq20112017(k,10);
        
        comptable(i,1)=date;
        comptable(i,2)=stn;
        comptable(i,3)=num2str(dis);
        comptable(i,4)=num2str(mag);
        comptable(i,5)=num2str(mag^3/dis);
        
        i=i+1;
        
    end
   
end

comptable = sortrows(comptable,-5);
