%% Purpose
% Search for precursors automatically for all detectable earthquakes with acceptable
% azimuthal angles accuracy by scanning 2 or 4 night time periods and all
% frequencies in ULF
% Divide ULF range to 10 equal ranges: 0.01 - 0.02 Hz, 0.02 - 0.03 Hz etc...
% Instead of using spectrogram and SSTF like Adib_16, this uses Hilbert
% transform and polarization ellipse
%% Loading universal variables
load VARIABLES_WORLD

mag_min=5.0;
dis_max=180;
depth_max=200;

LT_start=[22,00,02,04];
LT_end=[00,02,04,06];
 
no_records=0;

%% Analysis starts
for  stn_num=1:size(stn_vec,1)
    for year=2005:2018
        %% Checking files availability
        
        stn=string(stn_vec(stn_num));
        filepathname=strcat('E:\Study\MAGDAS DATA\',stn,'\',stn,string(year),'S.mat');
        filecheck=exist(filepathname,'file');
        
        if filecheck~=2 continue; end
        
        disp(sprintf('Running on %s%d..',stn,year));
        
        datenum_start=datenum([year,01,01,00,00,00]);
        datenum_end=datenum([year,12,31,23,59,59]);
        
        stn_latlon=[stn_MAGDAS(stn_num,2:3)];
        region=string(station_vec(stn_num,2));
        %% Building earthquakes table
        EQ_table=NaN(1,11);
        j=1;
        for i=1:length(EQ_WORLD)
            if datenum(EQ_WORLD(i,1:6))>=datenum_start && datenum(EQ_WORLD(i,1:6))<=datenum_end
                EQ_table(j,:)=EQ_WORLD(i,:);
                j=j+1;
            end
            if datenum(EQ_WORLD(i,1:6))>datenum_end
                break;
            end
        end
        
        i=1;
        EQ_sel=NaN(1,8);
        for j=1:size(EQ_table,1)
            if  EQ_table(j,10)>=mag_min
                EQ_latlon=[EQ_table(j,7:8)];
                EQ_dis=deg2km(distance('gc',stn_latlon,EQ_latlon));
                EQ_time=datenum(EQ_table(j,1:6));
                EQ_mag=EQ_table(j,10);
                EQ_Ks=(10^(0.75*EQ_mag))/(EQ_dis+100);
                EQ_depth=EQ_table(j,9);
                EQ_angle=azimuth('gc',stn_latlon(1:2),EQ_latlon(1:2),'degree');
                
                if EQ_dis<=dis_max && EQ_depth<=depth_max
                    EQ_sel(i,1)=EQ_time;
                    EQ_sel(i,2)=EQ_dis;
                    EQ_sel(i,3)=EQ_mag;
                    EQ_sel(i,4)=EQ_Ks;
                    EQ_sel(i,5)=EQ_depth;
                    EQ_sel(i,6:7)=EQ_latlon;
                    EQ_sel(i,8)=EQ_angle;
                    i=i+1;
                end
            end
        end
        
        EQ_selselidx=find( 0.025*EQ_sel(:,2)<EQ_sel(:,3)-4.5 | EQ_sel(:,3)>=6.5 | EQ_sel(:,2)<=30 );
        EQ_selsel=EQ_sel(EQ_selselidx,:);
        
        if isempty(EQ_selselidx)
            disp('No detectable/significant earthquakes. Skipping..');
            pause(0.5);
            clc;
            continue;
        else
            disp(sprintf('Found %d earthquake(s) to be analyzed',numel(EQ_selselidx)));
        end
        %% Loading geomagnetic data files
        days_num=round(datenum_end-datenum_start);
        if EQ_selsel(1,1)-datenum_start<60
            disp('Previous year file need to be loaded..')
            filepathname2=strcat('E:\Study\MAGDAS DATA\',stn,'\',stn,string(year-1),'S.mat');
            if exist(filepathname2,'file')
                matname(1)=strcat(stn,string(year-1),'S');
                matname(2)=strcat(stn,string(year),'S');
                A=load(matname(1));
                B=load(matname(2));
                H=vertcat(A.H,B.H);
                D=vertcat(A.D,B.D);
                Z=vertcat(A.Z,B.Z);
                UT1m=horzcat(A.UT1m,B.UT1m);
                
                H=H((end-(days_num+61)*86400)+1:end);
                D=D((end-(days_num+61)*86400)+1:end);
                Z=Z((end-(days_num+61)*86400)+1:end);
                UT1m=UT1m((end-(days_num+61)*86400)+1:end);
                
                disp('Previous year file has been loaded.');
                prevyear=1;
            else
                matname=strcat(stn,string(year),'S');
                load(matname);
                H_nan=NaN(61*86400,1); D_nan=NaN(61*86400,1);
                Z_nan=NaN(61*86400,1); UT1m_nan=UT1m(1)-61:1/86400:UT1m(1)-(1/86400);
                H=vertcat(H_nan,H); D=vertcat(D_nan,D);
                Z=vertcat(Z_nan,Z); UT1m=horzcat(UT1m_nan,UT1m);
                
                disp('Previous year file is unavailable. Proceeding..');
                prevyear=2;
            end
            
        else
            matname=strcat(stn,string(year),'S');
            load(matname);
            prevyear=0;
        end
        
        datadatenum_start=min(UT1m); datadatenum_end=max(UT1m);
        datenum_vec=datadatenum_start:floor(datadatenum_end);
        days_vec=datetime(datenum_vec,'Format','dd/MM/yyyy','ConvertFrom','datenum');
        days_num=length(days_vec);
        
        date_vec=datevec(datenum_vec);
        year_vec=unique(date_vec(:,1));
        if prevyear==2 year_vec=year_vec(end); end
        %% Removing noise/outlier
        H_dn=H; D_dn=D; Z_dn=Z;
        G_dn=sqrt(H_dn.^2+D_dn.^2);
        
        for i=1:5
            H_dn=medfilt1(H_dn,30,'omitnan','truncate');
            D_dn=medfilt1(D_dn,30,'omitnan','truncate');
            Z_dn=medfilt1(Z_dn,30,'omitnan','truncate');
            G_dn=medfilt1(G_dn,30,'omitnan','truncate');
        end
        
        for i=1:5
            
            H_sig=std(H_dn,'omitnan');
            H_mu=mean(H_dn,'omitnan');
            D_sig=std(D_dn,'omitnan');
            D_mu=mean(D_dn,'omitnan');
            Z_sig=std(Z_dn,'omitnan');
            Z_mu=mean(Z_dn,'omitnan');
            G_sig=std(G_dn,'omitnan');
            G_mu=mean(G_dn,'omitnan');
            
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
                if G_dn(j)>G_mu+5*G_sig||G_dn(j)<G_mu-5*G_sig
                    G_dn(j)=NaN;
                end
            end
        end
        %% Nighttime data in LT
        for l=1:numel(LT_start)
            LT_length(l)=LT_end(l)-LT_start(l);
            if LT_length(l)<0 LT_length(l)=LT_length(l)+24; end
            t_zone=timezone(stn_latlon(2),'degrees');
            t_start(l)=LT_start(l)+t_zone;
            if t_start(l)>24 t_start(l)=t_start(l)-24; end
            if t_start(l)<0 t_start(l)=t_start(l)+24;end
            t_starts(l)=t_start(l)*3600;
            t_int(l)=LT_length(l)*3600;
            
            for i=1:days_num
                H_night(:,i,l)=H_dn(t_starts(l)+1:t_starts(l)+t_int(l));
                D_night(:,i,l)=D_dn(t_starts(l)+1:t_starts(l)+t_int(l));
                Z_night(:,i,l)=Z_dn(t_starts(l)+1:t_starts(l)+t_int(l));
                G_night(:,i,l)=G_dn(t_starts(l)+1:t_starts(l)+t_int(l));
                UT1m_night(:,i,l)=UT1m(t_starts(l)+1:t_starts(l)+t_int(l));
                t_starts(l)=t_starts(l)+86400;
            end
        end 
        %% Spectrogram calculation
        n_win=1800;
        n_ovrlp=0.5*n_win;
        n_fft=n_win;
        fs=1;
        
        for l=1:numel(LT_start)
            for i=1:days_num
                [Z_spec,F_spec]=pwelch(Z_night(:,i,l),hamming(n_win),n_ovrlp,n_fft,fs);
                [G_spec,F_spec]=pwelch(G_night(:,i,l),hamming(n_win),n_ovrlp,n_fft,fs);
                
                bin001=find(round(F_spec,4)==0.0100); bin01=find(round(F_spec,4)==0.1000);
                
                Z_spec1=Z_spec(bin001:bin01,:);
                G_spec1=G_spec(bin001:bin01,:);
                F_spec1=F_spec(bin001:bin01);
                
                f_rangeint=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10];
                for a=1:numel(f_rangeint) f_rangeidx(a)=find(F_spec1==f_rangeint(a)); end
                
                for a=1:numel(f_rangeint)-1
                    for b=1:size(G_spec1,2)
                        Z_spec_band(:,a)=Z_spec1(f_rangeidx(a):f_rangeidx(a+1)-1);
                        G_spec_band(:,a)=G_spec1(f_rangeidx(a):f_rangeidx(a+1)-1);
                    end
                    F_spec_band(:,a)=F_spec1(f_rangeidx(a):f_rangeidx(a+1));
                end
            
                ZG(:,i,l)=mean(Z_spec_band./G_spec_band,1);
                
                H_pass=filter1('bp',H_night(:,i,l),'fc',[0.01 0.1],'fs',fs);
                D_pass=filter1('bp',D_night(:,i,l),'fc',[0.01 0.1],'fs',fs);
                
                H_hilb=hilbert(H_pass); D_hilb=hilbert(D_pass);
                
                PE_top=2*abs(H_hilb).*abs(D_hilb).*cos(angle(H_hilb)-angle(D_hilb));
                PE_bot=(abs(H_hilb)).^2-(abs(D_hilb)).^2;
                
                PE_theta=(atan2(PE_top,PE_bot))/2;
                PE_alpha1=pi-PE_theta;
                PE_alpha2=ones(1,length(PE_alpha1));
                for a=1:numel(PE_alpha1)
                   if PE_alpha1(a)==pi/2
                       PE_alpha2(a)=3*pi/2;
                   elseif PE_alpha1(a)>pi/2 && PE_alpha1(a)<pi
                       PE_alpha2(a)=PE_alpha1(a)+pi;
                   elseif PE_alpha1(a)==pi
                       PE_alpha2(a)=0;
                   elseif PE_alpha1(a)>pi && PE_alpha1(a)<3*pi/2
                       PE_alpha2(a)=PE_alpha1(a)-pi;
                   elseif PE_alpha1(a)==3*pi/2
                       PE_alpha2(a)=pi/2;
                   end
                end
                
                
                
            end
        end
        %% Identifying disturbed days
        geomag_day1=find(datenum(index_geomag(:,1),1,index_geomag(:,2))==datadatenum_start);
        for l=1:numel(LT_start)
            geomag_day=geomag_day1(t_start(l));    %+1 hour backward
            k=1;
            for i=geomag_day:24:length(index_geomag)
                ap(:,k,l)=index_geomag(i:i+3,6);   %+1 hour onward
                Dst(:,k,l)=index_geomag(i:i+3,5);
                k=k+1;
                if k>days_num break; end
            end
        end
        
        ZG_quiet=ZG;
        for l=1:numel(LT_start)
            for i=1:days_num
                if any(ap(:,i,l)>50) || any(Dst(:,i,l)<-50)
                    ZG_quiet(:,i,l)=deal(NaN);
                end
            end
        end
        
        if all(all(isnan(ZG_quiet))) 
            disp('No data available prior to the earthquake(s). Skipping');
            pause(0.5);
            continue; 
        end
        %% Identifying precursors
        for a=1:size(EQ_selsel,1)
            for b=1:24
                EQ_foundpre{a,b}='N/A';
            end
        end
        
        days_backwards=[30,45,60];
        for n=1:size(EQ_selsel,1)
            
            time_now=datetime('now','Format','dd-MM-yyyy HH:mm');
            if 0.025*EQ_selsel(n,2)<EQ_selsel(n,3)-4.5 EQ_dtctblty='Yes';
            else EQ_dtctblty='No'; end
            
            EQ_foundpre{n,1}=char(time_now);
            EQ_foundpre{n,2}=sprintf('%s: %d %d %d)',stn,year_vec);
            EQ_foundpre{n,3}=sprintf('%s',region);
            EQ_foundpre{n,4}=[];
            EQ_foundpre{n,5}=sprintf('%s',datetime(EQ_selsel(n,1),'ConvertFrom','datenum','Format','dd-MM-yyyy HH:mm:ss'));
            EQ_foundpre{n,6}=sprintf('%.1f',EQ_selsel(n,3));
            EQ_foundpre{n,7}=sprintf('%.0f',EQ_selsel(n,5));
            EQ_foundpre{n,8}=sprintf('%.0f',EQ_selsel(n,2));
            EQ_foundpre{n,9}=sprintf('%.2f',EQ_selsel(n,4));
            EQ_foundpre{n,10}=sprintf('%s',EQ_dtctblty);
            EQ_foundpre{n,11}=sprintf('%.2f',EQ_selsel(n,8));
            EQ_foundpre{n,12}=sprintf('%.2f',EQ_selsel(n,6));
            EQ_foundpre{n,13}=sprintf('%.2f',EQ_selsel(n,7));
            EQ_foundpre{n,14}=[];
            
            for f=1:numel(days_backwards)
                
                days_bckwrd=days_backwards(f);
                
                ZG_flex=ZG_quiet;
                EQ_dateidx=find(floor(EQ_selsel(n,1))==datenum_vec);
                foundpre=0;
                days_range=EQ_dateidx:-1:EQ_dateidx-days_bckwrd;
                
                for m=1:10
                    
                    p=0;
                    clearvars max_ZG max_ZGidx
                    for o=EQ_dateidx:-1:EQ_dateidx-days_bckwrd
                        p=p+1;
                        max_ZG(p)=max(max(ZG_flex(:,o,:)));
                        max_ZGlinidx=find(ZG_flex(:,o,:)==max_ZG(p));
                        mat_dim=[size(ZG_flex,1),size(ZG_flex,3)];
                        if ~isempty(max_ZGlinidx)
                            max_ZGlinidx=max_ZGlinidx(1);
                            [max_ZGidx(1,p),max_ZGidx(2,p)]=ind2sub(mat_dim,max_ZGlinidx);
                        else
                            max_ZGidx(1:2,p)=[NaN,NaN]; end
                    end
                    
                    max_ZG_flex=max_ZG;
                    for q=1:size(max_ZG,2)
                        
                        [r,s]=max(max_ZG_flex);
                        
                        if isnan(r) continue; end
                        
                        azim_ZG=azim_theta(max_ZGidx(1,s),days_range(s),max_ZGidx(2,s));
                        azim_EQ=EQ_selsel(n,8);
                        
                        diff_angle=azim_ZG-azim_EQ;
                        if diff_angle<-180 diff_angle=diff_angle+360;
                        elseif diff_angle>180 diff_angle=diff_angle-360; end
                        if abs(diff_angle)<=20
                            EQ_foundpre{n,15}=sprintf('%s',days_vec(days_range(s)));
                            EQ_foundpre{n,16}=sprintf('%02d - %02d',LT_start(max_ZGidx(2,s)),LT_end(max_ZGidx(2,s)));
                            EQ_foundpre{n,17}=sprintf('%.4f',r);
                            EQ_foundpre{n,18}=sprintf('%.4f',Z_spec_mu(max_ZGidx(1,s),days_range(s),max_ZGidx(2,s)));
                            EQ_foundpre{n,19}=sprintf('%.4f',G_spec_mu(max_ZGidx(1,s),days_range(s),max_ZGidx(2,s)));
                            EQ_foundpre{n,20}=sprintf('%d',floor(EQ_selsel(n,1))-datenum_vec(days_range(s)));
                            EQ_foundpre{n,21}=sprintf('%.4f',round(F_spec_mu1(max_ZGidx(1,s)),4));
                            EQ_foundpre{n,22}=sprintf('%.2f',azim_ZG);
                            EQ_foundpre{n,23}=sprintf('%.2f',diff_angle);
                            EQ_foundpre{n,24}=sprintf('%d',m);
                            foundpre=1;
                            break;
                        else
                            max_ZG_flex(s)=NaN;
                            if all(isnan(max_ZG_flex))
                                break;
                            end
                        end
                    end
                    
                    if foundpre==1
                        break;
                    else
                        for x=1:length(days_range)
                            if isnan(max_ZGidx(1,x)) continue; end
                            ZG_flex(max_ZGidx(1,x),days_range(x),max_ZGidx(2,x))=NaN;
                        end
                    end
                end
                if foundpre==1 break; end
            end
        end
        %% Writing to xls file
        xlsfilename='D:\Study\Record of appearance of earthquakes geomagnetic precursor.xlsx';
        xlssheetname='Adib_14g';
        
        [readxcel1,readxcel2,readxcel]=xlsread(xlsfilename,xlssheetname);
        empcell=size(readxcel,1)+1;
        xlsrange=strcat('A',string(empcell),':X',string(empcell+size(EQ_foundpre,1)-1));
        xlswrite(xlsfilename,EQ_foundpre,xlssheetname,xlsrange)
        
        disp(sprintf('Records made. Analysis on %s%d has completed.',stn,year));
        %% Clearing variables
        no_records=no_records+size(EQ_foundpre,1);
        clearvars -except mag_min dis_max depth_max LT_start LT_end stn_num year stn_vec stn_MAGDAS station_vec index_geomag EQ_WORLD no_records
        pause(0.5);
        clc;
    end
end

disp(sprintf('Analysis ends. %d of earthquakes recorded.',no_records));
