%% Purpose
%To test Adib_15 and runs automatically
%% Read xls file and load universal variables
xlsfilename='D:\Study\Record of appearance of earthquakes geomagnetic precursor.xlsx';
xlssheetname='Adib_15';
[readxcel1,readxcel2,readxcel]=xlsread(xlsfilename,xlssheetname);
today=strcat(char(datetime('today','Format','dd-MM-yyyy')),'\Adib_15B_AUTO\');
load VARIABLES_WORLD
%% Start of loop
for precursor=2:size(readxcel,1)
     
    if strcmp(readxcel{precursor,15},'N/A') continue; end
    
    %% Place and time
    stn1=readxcel{precursor,2};
    stn=stn1(1:3);
    date_EQ1=readxcel{precursor,5};
    date_EQ=date_EQ1(1:9);
    datenum_EQ=datenum(date_EQ); 
    
    datenum_start60=datenum_EQ-60;
    datenum_end30=datenum_EQ+30;
    
    year_array=unique(datetime([datenum_start60,datenum_end30],'ConvertFrom','datenum','Format','yyyy'));
    filepathname=strcat('E:\Study\MAGDAS DATA\',stn,'\',stn,string(year_array(:)),'S.mat');
  
    for i=1:size(year_array,2)
        filecheck=exist(filepathname(i),'file');
        if (i==1 && filecheck==2) || numel(year_array)==1
            date_start=datevec(datenum_start60);
        elseif i==1 && filecheck~=2
            date_start=[str2num(sprintf('%s',year_array(i+1))),01,01];
        end
        if (i==2 && filecheck==2) || numel(year_array)==1
            date_end=datevec(datenum_end30);;
        elseif i==2 && filecheck~=2
            date_end=[str2num(sprintf('%s',year_array(i-1))),12,31];
        end
    end

    disp(sprintf('Analyzing data of %s from %s to %s..',stn,datetime(date_start,'Format','dd/MM/yyyy'),datetime(date_end,'Format','dd/MM/yyyy')));
    %% Customization
    mag_min=5.0;                     
    dis_max=180;                   
    depth_max=200;                   
    %% Frequency band and LT night time
    %Frequency band
    f_1=readxcel{precursor,21};
    f_2=f_1;
    f_array=unique([f_1,f_2]);
    
    %Nighttime setup
    LT_read=readxcel{precursor,16};
    LT_start=str2num(LT_read(1:2));
    LT_end=str2num(LT_read(end-1:end));
    %% Time period
    date_start1=datenum(horzcat(date_start(1:3),[0,0,0]));
    date_end1=datenum(horzcat(date_end(1:3),[23,59,59]));
    datenum_start=datenum(date_start);
    datenum_end=datenum(date_end);
    days_num=datenum_end-datenum_start+1;
    days_vec=datetime(datevec(datenum_start:1:datenum_end),'Format','dd/MM/yyyy');
    
    if date_start(1)==date_end(1)
        year=string(date_start(1));
    else
        year=sprintf('%s-%s',string(date_start(1)),string(date_end(1)));
    end
    
    j=date_start(1);
    for i=1:3
        year_vec(1,i)=j;
        if j==date_end(1)
            break
        end
        j=j+1;
    end
    %% Loading files
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
    %% Station setting
    for i=1:length(stn_MAGDAS)
        if strcmp(stn,stn_vec(i))
            stn_num=i;
            station=string(station_vec(i,1));
            region=string(station_vec(i,2));
            break;
        end
    end
    
    stn_latlon=[stn_MAGDAS(stn_num,2:3)];
    %% Building earthquakes table
    EQ_table=NaN(1,11);
    j=1;
    for i=1:length(EQ_WORLD)
        if datenum(EQ_WORLD(i,1:6))>=datenum_start && datenum(EQ_WORLD(i,1:6))<=datenum_end
            EQ_table(j,:)=EQ_WORLD(i,:);
            j=j+1;
        end
        if datenum(EQ_WORLD(i,1:6))>datenum_end
            break
        end
    end
    
    i=1;
    EQ_sel=NaN(1,11);
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
    colpal=jet(size(EQ_sel,1));
    EQ_sel(:,9:11)=colpal;
    %% Geomagnetic data pre-processing
    data_start=round(min(UT1m));
    data_end=floor(max(UT1m));
    
    day_start=datenum_start-data_start+1;
    day_end=floor(datenum_end)-data_start+1;
    datenum_vec=datenum_start:datenum_end;
    UT1m_sel=datevec(UT1m(86400*day_start-86399:86400*day_end));
    
    H_dn=H(86400*day_start-86399:86400*day_end);
    D_dn=D(86400*day_start-86399:86400*day_end);
    Z_dn=Z(86400*day_start-86399:86400*day_end);
    G_dn=sqrt(H_dn.^2+D_dn.^2);
    
    %Removing noise/outlier
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
        H_night(:,i)=H_dn(t_starts+1:t_starts+t_int);
        D_night(:,i)=D_dn(t_starts+1:t_starts+t_int);
        Z_night(:,i)=Z_dn(t_starts+1:t_starts+t_int);
        G_night(:,i)=G_dn(t_starts+1:t_starts+t_int);
        t_starts=t_starts+86400;
    end
    %% ap and Dst indices
    geomag_day1=find(datenum(index_geomag(:,1),1,index_geomag(:,2))==datenum_start);
    geomag_day=geomag_day1(t_start);    %+1 hour backward
    k=1;
    for i=geomag_day:24:length(index_geomag)
        ap(:,k)=index_geomag(i:i+3,6);   %+1 hour onward
        Dst(:,k)=index_geomag(i:i+3,5);
        k=k+1;
        if k>days_num break; end
    end
    %% PSD calculation
    n_win=900;
    n_ovrlp=[];
    n_fft=n_win;
    fs=1;
    
    for i=1:days_num
    
        H_nightseg=reshape(H_night(:,i),n_win,[]);
        D_nightseg=reshape(D_night(:,i),n_win,[]);
        Z_nightseg=reshape(Z_night(:,i),n_win,[]);
        G_nightseg=reshape(G_night(:,i),n_win,[]);
        
        for m=1:size(H_nightseg,2)
            [H_spec(:,m),F_spec]=pwelch(H_nightseg(:,m),hamming(n_win),n_ovrlp,n_fft,fs);
            [D_spec(:,m),F_spec]=pwelch(D_nightseg(:,m),hamming(n_win),n_ovrlp,n_fft,fs);
            [Z_spec(:,m),F_spec]=pwelch(Z_nightseg(:,m),hamming(n_win),n_ovrlp,n_fft,fs);
            [G_spec(:,m),F_spec]=pwelch(G_nightseg(:,m),hamming(n_win),n_ovrlp,n_fft,fs);
        end
                
        Z_spec_mu=mean(abs(Z_spec),2);
        G_spec_mu=mean(abs(G_spec),2);
        
        bin_1=find(round(F_spec,4)==round(f_1,4));
        bin_2=find(round(F_spec,4)==round(f_2,4));
        
        Z_spec_f=mean(Z_spec_mu(bin_1:bin_2));
        G_spec_f=mean(G_spec_mu(bin_1:bin_2));
        
        ZG(i)=Z_spec_f/G_spec_f;
        
        XX=real(H_spec(bin_1:bin_2,:));
        YY=real(D_spec(bin_1:bin_2,:));
        ZZ=real(Z_spec(bin_1:bin_2,:));
        XY_mat=[XX;YY];
        ZZ_mat=(ZZ');
        AB=(inv(XY_mat*XY_mat'))*(XY_mat*ZZ_mat);
        azim_amp(i)=sqrt(AB(1)^2+AB(2)^2);
        azim_theta1=atan2d(AB(2),AB(1));
        if azim_theta1<0 azim_theta1=360+azim_theta1; end
        azim_theta(i)=azim_theta1;
    end
    %% Disturbed days and anomalies detection
    %mu+/-2sigma
    ZG_mu=nanmean(ZG);
    ZG_sig=nanstd(ZG);
    ZG_mup2sig=ones(1,days_num).*(ZG_mu+2*ZG_sig);
    % ZG_mum2sig=ones(1,days_num).*(ZG_mu-2*ZG_sig);
    
    dstb=NaN(1,days_num);
    anom=NaN(1,days_num);
    
    for i=1:days_num
        if ZG(i)>ZG_mup2sig(1) %|| ZG(i)<ZG_mum2sig(1)
            if  any(ap(:,i)>50) || any(Dst(:,i)<-50)
                dstb(1,i)=ZG(i);
            else
                anom(1,i)=ZG(i);
            end
        end
    end
    %% Specifying date of interest
    DOI_0=datenum(date_EQ1);
    DOI_dates=datetime([DOI_0,DOI_0-3.0,DOI_0-7.0,DOI_0-14.0,DOI_0-30.0],'ConvertFrom','datenum','Format','dd/MM/yyyy HH');
    DOI_labels={'3D','1W','2W','1M'};
    %% Figure 1
    if numel(f_array)>1
        title_sup=sprintf('%s station, %s - %s | %cf = %.4f - %.4f Hz, %cLT = %02d - %02d hr',station,datetime(datevec(datenum(date_start)),'Format','dd/MM/yyyy'),datetime(datevec(datenum(date_end)),'Format','dd/MM/yyyy'),char(916),f_1,f_2,char(916),LT_start,LT_end);
    else
        title_sup=sprintf('%s station, %s - %s | f = %.4f Hz, %cLT = %02d - %02d hr',station,datetime(datevec(datenum(date_start)),'Format','dd/MM/yyyy'),datetime(datevec(datenum(date_end)),'Format','dd/MM/yyyy'),f_1,char(916),LT_start,LT_end); 
    end
    
    f1=figure(1);
    stitle=title(sprintf('%s',title_sup));
    stitle_pos =get(stitle,'position');
    stitle_pos=[stitle_pos(1) stitle_pos(2)+0.03 stitle_pos(3)];
    set(stitle,'position',stitle_pos);
    set(gca,'XTickLabel',' ');
    
    tightaxes=tight_subplot(2,1,0.005,[0.1 0.05]);     % (Nh, Nw, gap, marg_h, marg_w)
    
    sp1=subplot(tightaxes(1));
    purply=[0.4940 0.1840 0.5560]; greeny=[0.4660 0.6740 0.1880];
    plot(days_vec,ap(1,:),'Color',purply);
    hold on
    plot(days_vec,Dst(1,:),'Color',greeny);
    plot(days_vec,ap(2:4,:),'Color',purply,'HandleVisibility','off');
    plot(days_vec,Dst(2:4,:),'Color',greeny,'HandleVisibility','off');
    plot(days_vec,ones(length(days_vec),1)*(50),'--','Color',purply,'HandleVisibility','off');
    plot(days_vec,ones(length(days_vec),1)*(-50),'--','Color',greeny,'HandleVisibility','off');
    hold off
    ylabel('Geomagnetic indices (nT)');
    legend({'ap','Dst'},'location','northwest','orientation','horizontal');
    if ~isnan(DOI_0)
        text(datetime(DOI_0,'ConvertFrom','datenum'),50,' Storm','Color',purply,'VerticalAlignment','top','HorizontalAlignment','left','FontSize',9);
        text(datetime(DOI_0,'ConvertFrom','datenum'),-50,' Storm','Color',greeny,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',9);
    end
    if max(max(ap))>50 && min(min(Dst))<-50 ylim([min(min(Dst))-10 max(max(ap))+10])
    elseif max(max(ap))>50 ylim([-60 max(max(ap))+10])
    elseif min(min(Dst))<-50 ylim([min(min(Dst))-10 60])
    else ylim([-60 60])
    end
    
    yyaxis right
    set(gca, 'YColor', 'k');
    if ~any(isnan(EQ_sel))
        for i=1:size(EQ_sel,1)
            hold on
            j=datetime(datevec(EQ_sel(i,1)),'Format','dd/MM/yyyy HH:mm:ss');
            label=sprintf('  %s, M%0.1f, D=\t%.0fKM, d=%.0fKM',string(datetime(j,'Format','dd/MM/yyyy')),EQ_sel(i,3),EQ_sel(i,5),EQ_sel(i,2));
            if 0.025*EQ_sel(i,2)<EQ_sel(i,3)-4.5
                plot(j,EQ_sel(i,4),'^','MarkerEdgeColor',EQ_sel(i,9:11),'MarkerFaceColor',EQ_sel(i,9:11),'MarkerSize',5,'HandleVisibility','off');
                text(j,EQ_sel(i,4),label,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',7.5);
            elseif EQ_sel(i,3)>6.5 && ~(0.025*EQ_sel(i,2)<EQ_sel(i,3)-4.5)
                plot(j,EQ_sel(i,4),'s','MarkerEdgeColor',EQ_sel(i,9:11),'MarkerFaceColor',EQ_sel(i,9:11),'MarkerSize',5,'HandleVisibility','off');
                text(j,EQ_sel(i,4),label,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',7.5);
            else
                plot(j,EQ_sel(i,4),'o','MarkerEdgeColor',EQ_sel(i,9:11),'MarkerFaceColor',EQ_sel(i,9:11),'MarkerSize',5,'HandleVisibility','off');
            end
            text(j,EQ_sel(i,4),sprintf('%.0f%c  ',EQ_sel(i,8),char(176)),'VerticalAlignment','middle','HorizontalAlignment','right','FontSize',7.5);
        end
        
        DOI_dates=datetime([DOI_0,DOI_0-3.0,DOI_0-7.0,DOI_0-14.0,DOI_0-30.0],'ConvertFrom','datenum','Format','dd/MM/yyyy HH');
        DOI_labels={'3D','1W','2W','1M'};
        plot([DOI_dates(1) DOI_dates(1)],[-1e5 1e5],'-','Color',[0 0 0],'HandleVisibility','off');
    end
    
    hold off
    set(gca,'XTickLabel','');
    xlim([min(days_vec) max(days_vec)]);
    if ~any(isnan(EQ_sel(1,:)))  ylim([0 1.1*max(EQ_sel(:,4))]); end
    ylabel('K_{LS} (a.u.)');
    
    
    sp2=subplot(tightaxes(2));
    bar(days_vec,ZG,'HandleVisibility','off');
    hold on
    plot(days_vec,ZG_mup2sig,'r--','HandleVisibility','off');
    if ~isnan(DOI_0)
        text(datetime(DOI_0,'ConvertFrom','datenum'),ZG_mup2sig(1),' \mu+2\sigma','Color','r','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',9);
    end
    % plot(days_vec,ZG_mum2sig,'k--','HandleVisibility','off');
    bar(days_vec,anom,'r');
    bar(days_vec,dstb,'g');
    if  ~any(isnan(EQ_sel(1,:))) 
        for i=1:5
            if i==1
                plot([DOI_dates(i) DOI_dates(i)],[-1e5 1e5],'-','Color',[0 0 0],'HandleVisibility','off');
            else
                plot([DOI_dates(i) DOI_dates(i)],[-1e5 1e5],':','Color',[0 0 0],'HandleVisibility','off');
                text(DOI_dates(i),1.2*max(ZG),DOI_labels(i-1),'VerticalAlignment','top','HorizontalAlignment','right','FontSize',7.5)
            end
        end
    end
    for v=1:days_num
        if ~isnan(anom(v))
            text(days_vec(v),ZG(v),sprintf('%.1f%c\n%.2f',azim_theta(v),char(176),ZG(v)),'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',7.5,'Color','r');
        end
    end
    hold off
    ylabel('Z/G (a.u.)');
    legend({'Anomalies','Disturbed'},'location','northwest','orientation','horizontal');
    xlim([min(days_vec) max(days_vec)]);
    ylim([0 1.2*max(ZG)])
    
    yyaxis right
    set(gca, 'YColor', 'k');
    plot(days_vec,azim_theta,'ko-.','HandleVisibility','off');
    xlim([min(days_vec) max(days_vec)]);
    ylim([0 360]);
    ylabel(sprintf('Azimuthal angle (%c)',char(176)));
    
    linkaxes([sp1,sp2],'x')
    
    set(gcf, 'Position', get(0, 'Screensize'));
    disp('Figure generated')
    %% Saving figures
    ddMMyy1=upper(string(datetime(min(days_vec),'Format','ddMMyy')));
    ddMMyy2=upper(string(datetime(max(days_vec),'Format','ddMMyy')));
    
    path1='D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\';
    mkdir(fullfile(path1,today));
    
    path2=strcat('D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\',today);
    
    figname1=char(strcat(stn,'_',ddMMyy1,'-',ddMMyy2));
    
    for i=1:100
        if exist(strcat(path2,figname1,sprintf('-%d',i),'.tif'))~=2
            figname1ext=strcat(figname1,char(sprintf('-%d',i)));
            saveas(f1,fullfile(path2,figname1ext),'tiff');
            figpath=strcat(path2,figname1ext,'.tif');
            break;
        end
    end
    %% Evaluate output
%     date_precursor=datetime(readxcel{precursor,15},'InputFormat','dd-MMM-yy','Format','dd/MM/yyyy');
%     index_precursor=find(days_vec==date_precursor);
%     
%     error_amp=(abs(ZG(index_precursor)-readxcel{precursor,17})/readxcel{precursor,17})*100;
%     
%     diff_angle=azim_theta(index_precursor)-readxcel{precursor,22};
%     if diff_angle<-180 diff_angle=diff_angle+360;
%     elseif diff_angle>180 diff_angle=diff_angle-360; end
%     error_azim=(abs(diff_angle)/360)*100;
    %% Write to xls file
    xlswrite(xlsfilename,cellstr(figpath),xlssheetname,strcat('Z',string(precursor)));
    xlswrite(xlsfilename,round(ZG_mu,4),xlssheetname,strcat('AA',string(precursor)));
    xlswrite(xlsfilename,round(ZG_mup2sig,4),xlssheetname,strcat('AB',string(precursor)));
    
    disp('Records made.');
    pause(0.5);
    %% Clearing variables
    clearvars -except xlsfilename xlssheetname readxcel EQ_WORLD index_geomag station_vec stn_MAGDAS stn_vec precursor today
    clc
    close all
end