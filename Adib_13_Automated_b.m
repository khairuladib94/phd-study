%% Purpose
%Automation - run on station & year included in stnyr_sel
%Can choose to combine with previous or following year

%To plot daily Z/G and geomagnetic indices temporal evolution with
%customizable segmentation period (hourly as default)
%Perform Welch in all four frequency channels at the same time
%Standardize Z and G individually and compute ZG only if Z or G >0.1 based
%on Ida et al. 2008
%Plot data in bar instead of line
%% Universal parameters
% stnyr_sel={'AMA2009','AMA2010','ANC2007','ANC2013','ANC2014','ASB2011','ASB2012','ASB2013','CDO2010','CDO2011','CDO2013','CDO2015','CEB2011','CEB2012','CEB2015','DAV2007','DAV2008','DAV2010','DAV2011','DAV2012','DAV2013','DAV2014','DAV2015','GSI2015','HLN2009','HLN2013','ICA2013','ICA2017','KTB2007','KTB2009','KTB2017','KUJ2014','KUJ2015','LAQ2009','LGZ2012','LWA2014','LWA2015','LWA2016','MCQ2008','MND2007','MND2008','MND2014','MUT2011','MUT2012','MUT2015','ONW2008','ONW2010','ONW2011','ONW2012','ONW2013','ONW2014','ONW2015','PTK2008','PTK2013','PTK2015','PTK2016','SCN2012','SCN2016','SCN2018','TNO2010','TNO2011','TNO2013','TNO2015'};
stnyr_sel={'ANC2007','ASB2013','CDO2015','CEB2012','HLN2009','HLN2013','ICA2013','KTB2007','KUJ2014','KUJ2015','LGZ2012','LWA2014','LWA2015','LWA2016','ONW2008','ONW2010','ONW2011','ONW2012','ONW2015','PTK2008','PTK2013','PTK2016','TNO2011','TNO2015'};
%% 
for U=1:length(stnyr_sel)
    disp(sprintf('Running on %s...',stnyr_sel{U}));
    %% Place and time inputs
    
    stn=stnyr_sel{U}(1:3);
    if numel(stnyr_sel{U})==7
        yr_sel1=str2num(stnyr_sel{U}(4:7));
        yr_sel2=yr_sel1;
    elseif numel(stnyr_sel{U})==11
        yr_sel1=str2num(stnyr_sel{U}(4:7));
        yr_sel2=str2num(stnyr_sel{U}(8:11));
    end
    
    date_start=[yr_sel1,01,01];           %Insert custom start and end dates
    date_end=  [yr_sel2,12,31];           %Period spanning through 3 consecutive years is the maximum
    %% Customization
    mag_min=5.0;                     %Minimum magnitude of earthquakes to be considered
    dis_max=180;                     %Maximum epicentral distance from the station
    
    %Frequency channels
    f_1=0.0067;
    f_2=0.0100;
    
    f_3=0.0100;
    f_4=0.0220;
    
    f_5=0.0220;
    f_6=0.0500;
    
    f_7=0.0500;
    f_8=0.1000;
    
    N_seg=3*3600;                     %Segmentation length
%   pfit_order=0;                   %Order of polynomial fit for ZG detrending
    movmean_val=86400/N_seg;        %Moving mean of ZG and indices
    
    %Nighttime indicator patch
    LT_start=22;
    LT_end=04;
    %% Time period
    date_start1=datenum(horzcat(date_start,[0,0,0]));
    date_end1=datenum(horzcat(date_end,[23,59,59]));
    days_vec=datetime(datevec(date_start1:1/(86400/N_seg):date_end1),'Format','dd/MM/yyyy HH');
    date_vec=datevec(date_start1:1/(86400/N_seg):date_end1);
    datenum_start=datenum(date_start1);
    datenum_end=datenum(date_end1);
    days_num=length(days_vec);
    
    if date_start(1)==date_end(1)
        year=string(date_start(1));
    else
        year=sprintf('%s-%s',string(date_start(1)),string(date_end(1)));
    end
    
    load VARIABLES_WORLD
    
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
    %% Station setting
    for i=1:length(stn_MAGDAS)
        if strcmp(stn,stn_vec(i))
            stn_num=i;
            station=string(station_vec(i));
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
    colpal(:,1)=0:10:180;
    colpal(:,2:4)=hot(size(colpal,1));
    for j=1:size(EQ_table,1)
        if  EQ_table(j,10)>=mag_min
            EQ_latlon=[EQ_table(j,7:8)];
            EQ_dis=deg2km(distance('gc',stn_latlon,EQ_latlon));
            EQ_time=datenum(EQ_table(j,1:6));
            EQ_mag=EQ_table(j,10);
            EQ_Ks=(10^(0.75*EQ_mag))/(EQ_dis+100);
            EQ_depth=EQ_table(j,9);
            EQ_angle=azimuth('gc',stn_latlon(1:2),EQ_latlon(1:2),'degree');
            if EQ_angle>180
                EQ_angle=EQ_angle-360;
            end
            EQ_discol=colpal(find(colpal(:,1)==round(EQ_dis,-1)),2:4);
            
            if EQ_dis<=dis_max
                EQ_sel(i,1)=EQ_time;
                EQ_sel(i,2)=EQ_dis;
                EQ_sel(i,3)=EQ_mag;
                EQ_sel(i,4)=EQ_Ks;
                EQ_sel(i,5)=EQ_depth;
                EQ_sel(i,6:7)=EQ_latlon;
                EQ_sel(i,8)=EQ_angle;
                EQ_sel(i,9:11)=EQ_discol;
                i=i+1;
            end
        end
    end
    %% Geomagnetic data pre-processing
    data_start=round(min(UT1m));
    data_end=floor(max(UT1m));
    
    day_start=datenum_start-data_start+1;
    day_end=floor(datenum_end)-data_start+1;
    datenum_vec=datenum_start:datenum_end;
    
    H_period=H(86400*day_start-86399:86400*day_end);
    D_period=D(86400*day_start-86399:86400*day_end);
    Z_period=Z(86400*day_start-86399:86400*day_end);
    G_period=sqrt(H_period.^2+D_period.^2);
    
    %Removing noise/outlier
    H_dn=medfilt1(H_period,'omitnan');
    D_dn=medfilt1(D_period,'omitnan');
    Z_dn=medfilt1(Z_period,'omitnan');
    G_dn=medfilt1(G_period,'omitnan');
    
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
            
            if H_dn(j)>H_mu+10*H_sig||H_dn(j)<H_mu-10*H_sig
                H_dn(j)=NaN;
            end
            if D_dn(j)>D_mu+10*D_sig||D_dn(j)<D_mu-10*D_sig
                D_dn(j)=NaN;
            end
            if Z_dn(j)>Z_mu+10*Z_sig||Z_dn(j)<Z_mu-10*Z_sig
                Z_dn(j)=NaN;
            end
            if G_dn(j)>G_mu+10*G_sig||G_dn(j)<G_mu-10*G_sig
                G_dn(j)=NaN;
            end
            
        end
        
    end
    
    H_seg=reshape(H_dn,N_seg,[]);
    D_seg=reshape(D_dn,N_seg,[]);
    Z_seg=reshape(Z_dn,N_seg,[]);
    G_seg=reshape(G_dn,N_seg,[]);
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
    t_end=t_start+LT_length;
    %% ap and Dst indices
    k=1;
    for i=1:length(index_geomag)
        m=datenum(index_geomag(i,1),1,index_geomag(i,2));
        if m>=datenum_start && m<=datenum_end
            AP(k,1)=index_geomag(i,6);
            DST(k,1)=index_geomag(i,5);
            k=k+1;
        end
    end
    
    idx_num=N_seg/3600;
    if idx_num>1
        ap=reshape(AP,idx_num,[]);
        Dst=reshape(DST,idx_num,[]);
        ap=mean(ap)';
        Dst=mean(Dst)';
    elseif idx_num==1
        ap=AP;
        Dst=DST;
    elseif idx_num<1
        idx_mult=days_num/length(AP);
        j=1;
        for i=1:length(AP)
            ap(j:j+idx_mult-1,1)=AP(i,1);
            Dst(j:j+idx_mult-1,1)=DST(i,1);
            j=j+idx_mult;
        end
    end
    %% PSD calculation
    n_win=1800;
    n_ovrlp=0.5*n_win;
    n_fft=n_win;
    fs=1;
    
    f_bot=[f_1,f_3,f_5,f_7];
    f_top=[f_2,f_4,f_6,f_8];
    
    for i=1:4
        bin_1=[];
        bin_2=[];
        
        [H_wlch,f_wlch]=pwelch(H_seg,hamming(n_win),n_ovrlp,n_fft,fs);
        [D_wlch,f_wlch]=pwelch(D_seg,hamming(n_win),n_ovrlp,n_fft,fs);
        [Z_wlch,f_wlch]=pwelch(Z_seg,hamming(n_win),n_ovrlp,n_fft,fs);
        [G_wlch,f_wlch]=pwelch(G_seg,hamming(n_win),n_ovrlp,n_fft,fs);
        
        f_wlch4=round(f_wlch,4);
        bin_1=find(f_wlch4==f_bot(i));
        bin_2=find(f_wlch4==f_top(i));
        if isempty(bin_1) || isempty(bin_2)
            f_wlch3=round(f_wlch,3);
            bin_1=floor(mean(find(f_wlch3==f_bot(i))));
            bin_2=floor(mean(find(f_wlch3==f_top(i))));
        end
        
        H_wlchf=H_wlch(bin_1:bin_2,:);
        H_wlchmu(:,:,i)=mean(H_wlchf,1,'omitnan');
        D_wlchf=D_wlch(bin_1:bin_2,:);
        D_wlchmu(:,:,i)=mean(D_wlchf,1,'omitnan');
        Z_wlchf=Z_wlch(bin_1:bin_2,:);
        Z_wlchmu(:,:,i)=mean(Z_wlchf,1,'omitnan');
        G_wlchf=G_wlch(bin_1:bin_2,:);
        G_wlchmu(:,:,i)=mean(G_wlchf,1,'omitnan');
    end
    
    %Normalization based on months
    for l=1:4
        for i=1:size(year_vec,2)
            same_yr=find(date_vec(:,1)==year_vec(i));
            mon_vec=unique(date_vec(same_yr,2),'stable');
            for j=1:size(mon_vec,1)
                same_mon=find(date_vec(:,2)==mon_vec(j));
                H_norm(1,same_mon,l)=normalize(H_wlchmu(1,same_mon,l));
                D_norm(1,same_mon,l)=normalize(D_wlchmu(1,same_mon,l));
                G_norm(1,same_mon,l)=normalize(G_wlchmu(1,same_mon,l));
                Z_norm(1,same_mon,l)=normalize(Z_wlchmu(1,same_mon,l));
            end
        end
    end
    
%     H_norm=normalize(H_wlchmu);
%     D_norm=normalize(D_wlchmu);
%     Z_norm=normalize(Z_wlchmu);
%     G_norm=normalize(G_wlchmu);
%     
    G_norm(G_norm==0)=NaN;
    %% Power ratio
    for j=1:4
        for i=1:days_num
%             if Z_norm(1,i,j)>0.1 || G_norm(1,i,j)>0.1
                ZG(1,i,j)=Z_norm(1,i,j)/G_norm(1,i,j);
%             else
%                 ZG(1,i,j)=NaN;
%             end
        end
    end
    ZG=movmean(ZG,movmean_val,'omitnan');
    
    %mu+/-2sigma
    ZG_mu=nanmean(ZG);
    ZG_sig=nanstd(ZG);
    ZG_mup2sig=ones(1,days_num).*(ZG_mu+2*ZG_sig);
    ZG_mum2sig=ones(1,days_num).*(ZG_mu-2*ZG_sig);
    %% Disturbed days and anomalies detection
    dstb=NaN(1,days_num,4);
    anom=NaN(1,days_num,4);
    for j=1:4
        for i=1:days_num
            if ZG(1,i,j)>ZG_mup2sig(1,i,j) || ZG(1,i,j)<ZG_mum2sig(1,i,j)
                if ap(i,1)>50 || Dst(i,1)<-50
                    dstb(1,i,j)=ZG(1,i,j);
                    if i<days_num
                        dstb(1,i+1,j)=ZG(1,i+1,j);
                    end
                    if i>1
                        dstb(1,i-1,j)=ZG(1,i-1,j);
                    end
                else
                    anom(1,i,j)=ZG(1,i,j);
                    if i<days_num
                        anom(1,i+1,j)=ZG(1,i+1,j);
                    end
                    if i>1
                        anom(1,i-1,j)=ZG(1,i-1,j);
                    end
                end
            end
        end
    end
    
    largeEQ=EQ_sel(find(EQ_sel(:,3)==max(EQ_sel(:,3))),1); 
    largeEQ_3D=datetime(largeEQ-3.0,'ConvertFrom','datenum','Format','dd/MM/yyyy HH');
    largeEQ_1W=datetime(largeEQ-7.0,'ConvertFrom','datenum','Format','dd/MM/yyyy HH');
    largeEQ_2W=datetime(largeEQ-14.0,'ConvertFrom','datenum','Format','dd/MM/yyyy HH');
    largeEQ_1M=datetime(largeEQ-30.0,'ConvertFrom','datenum','Format','dd/MM/yyyy HH');
    %% Azimuthal angle calculation
    % XX=H_seg;
    % YY=D_seg;
    % ZZ=Z_seg;
    % n_win_dir=120;
    % n_ovrlp_dir=60;
    % n_fft_dir=n_win_dir;
    % for i=1:size(XX,2)
    %     XX_spec=spectrogram(XX(:,i),'yaxis',hamming(n_win_dir),n_ovrlp_dir,n_fft_dir,1);
    %     YY_spec=spectrogram(YY(:,i),'yaxis',hamming(n_win_dir),n_ovrlp_dir,n_fft_dir,1);
    %     ZZ_spec=spectrogram(ZZ(:,i),'yaxis',hamming(n_win_dir),n_ovrlp_dir,n_fft_dir,1);
    %     bin036=round(0.36/(0.5/size(XX_spec,1)));
    %     XX_036=real(XX_spec(bin036,:));
    %     YY_036=real(YY_spec(bin036,:));
    %     ZZ_036=real(ZZ_spec(bin036,:));
    %     XY_mat=[XX_036;YY_036];
    %     ZZ_mat=(ZZ_036');
    %     AB=real((inv(XY_mat*XY_mat'))*(XY_mat*ZZ_mat));
    %     AB_amp(1,i)=sqrt(AB(1)^2+AB(2)^2);
    %     AB_theta(1,i)=atan2d(AB(2),AB(1));
    % end
    %     XX=H_wlchf;
    %     YY=D_wlchf;
    %     ZZ=Z_wlchf;
    %     for i=1:size(XX,2)
    %         XY_mat=[XX(:,i) YY(:,i)]';
    %         ZZ_mat=ZZ(:,i);
    %         AB=real((inv(XY_mat*XY_mat'))*(XY_mat*ZZ_mat));
    %         AB_amp(1,i)=(sqrt(AB(1)^2+AB(2)^2));
    %         AB_theta(1,i)=atan2d(AB(2),AB(1));
    %     end
    %% Generate figure
    title_sup=sprintf('%s station, %s - %s',station,datetime(datevec(datenum(date_start)),'Format','dd/MM/yyyy'),datetime(datevec(datenum(date_end)),'Format','dd/MM/yyyy'));
    
    f1=figure(1);
    stitle=suptitle(sprintf('%s',title_sup));
    stitle_pos =get(stitle,'position');
    stitle_pos=[stitle_pos(1) stitle_pos(2)+0.01 stitle_pos(3)];
    set(stitle,'position',stitle_pos);
    set(gca,'XTickLabel','','YTickLabel','');
    
    tightaxes=tight_subplot(5,1,0,[0.1 0.05]);     % (Nh, Nw, gap, marg_h, marg_w)
    
    sp1=subplot(tightaxes(1));
    plot(days_vec,ap,'k');
    hold on
    plot(days_vec,Dst,'r');
    plot(days_vec,ones(length(days_vec),1)*(50),'k--');
    plot(days_vec,ones(length(days_vec),1)*(-50),'r--');
    hold off
    ylabel('Dst               ap');
    if max(ap)>50 && min(Dst)<-50 ylim([min(Dst)-10 max(ap)+10])
    elseif max(ap)>50 ylim([-60 max(ap)+10])
    elseif min(Dst)<-50 ylim([min(Dst)-10 60])
    else ylim([-60 60])
    end
    yyaxis right
    for i=1:size(EQ_sel,1)
        hold on
        j=datetime(datevec(EQ_sel(i,1)),'Format','dd/MM/yyyy HH:mm:ss');
        label={sprintf('%s',string(datetime(j,'Format','dd/MM/yyyy')));
                sprintf('M%0.1f D=\t%.0fKM',EQ_sel(i,3),EQ_sel(i,5));
                sprintf('d=%.0fKM',EQ_sel(i,2))};
        if 0.025*EQ_sel(i,2)<EQ_sel(i,3)-4.5
            plot(j,EQ_sel(i,3),'^','MarkerEdgeColor',EQ_sel(i,9:11),'MarkerFaceColor',EQ_sel(i,9:11),'MarkerSize',5);
            text(j,EQ_sel(i,3),label,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',9);
        else
            plot(j,EQ_sel(i,3),'o','MarkerEdgeColor',EQ_sel(i,9:11),'MarkerFaceColor',EQ_sel(i,9:11),'MarkerSize',5);
        end
        if EQ_sel(i,3)==max(EQ_sel(:,3))
            t=i;
            u=j;
            plot([u u],[-1e5 1e5],'k-');
        end
    end
    hold off
    set(gca,'YColor', 'b','XTickLabel','')
    xlim([min(days_vec) max(days_vec)])
    ylim([5.0 9.5])
    ylabel('Magnitude');
    
    
    sp_list=char({'sp2','sp3','sp4','sp5'});
    for i=1:4
        eval([sp_list(i,:) '=subplot(tightaxes(i+1));']);
        bar(days_vec,ZG(:,:,i))
        hold on
        for k=1:days_num/(86400/N_seg)
            xpatch=[k k+1 k+1 k];
            ypatch=[-1e+5 -1e+5 1e+5 1e+5];
            if isnan(Z_norm(:,k*(86400/N_seg),i)) || isnan(G_norm(:,k*(86400/N_seg),i)) 
                patch(xpatch,ypatch,'y','FaceAlpha',0.2,'EdgeAlpha',0);
            end
        end
        bar(days_vec,ZG(:,:,i),'b')
        bar(days_vec,anom(:,:,i),'r');
        bar(days_vec,dstb(:,:,i),'g');
        plot(days_vec,ZG_mum2sig(:,:,i),'k--',days_vec,ZG_mup2sig(:,:,i),'k--');
        legend({sprintf('\\Deltaf=%.3f-%.3f Hz',f_bot(i),f_top(i))},'location','northwest','orientation','horizontal');
        plot([u u],[-1e5 1e5],'k','HandleVisibility','off');
        plot([largeEQ_3D largeEQ_3D],[-1e5 1e5],'--','Color',[0.5 0.5 0.5],'HandleVisibility','off');
        plot([largeEQ_1W largeEQ_1W],[-1e5 1e5],'--','Color',[0.5 0.5 0.5],'HandleVisibility','off');
        plot([largeEQ_2W largeEQ_2W],[-1e5 1e5],'--','Color',[0.5 0.5 0.5],'HandleVisibility','off');
        plot([largeEQ_1M largeEQ_1M],[-1e5 1e5],'--','Color',[0.5 0.5 0.5],'HandleVisibility','off');

        if i==4
            text(largeEQ_3D,0,'-3D','VerticalAlignment','top','HorizontalAlignment','right','FontSize',9)
            text(largeEQ_1W,0,'-1W','VerticalAlignment','top','HorizontalAlignment','right','FontSize',9)
            text(largeEQ_2W,0,'-2W','VerticalAlignment','top','HorizontalAlignment','right','FontSize',9)
            text(largeEQ_1M,0,'-1M','VerticalAlignment','top','HorizontalAlignment','right','FontSize',9)
        end
        hold off
        ylabel('Z/G power ratio');
        if 0.8*min(ZG(:,:,i))<0 ylim([0.8*min(ZG(:,:,i)) 1.2*max(ZG(:,:,i))]);
        else ylim([-10 1.2*max(ZG(:,:,i))]); end
        xlim([min(days_vec) max(days_vec)])
        set(gca,'YTickLabelMode','auto')
        if i~=4 set(gca,'XTickLabel','');  end
    end
    
    linkaxes([sp1,sp2,sp3,sp4,sp5],'x')
    linkaxes([sp2,sp3,sp4,sp5],'y')
    
    set(gcf, 'Position', get(0, 'Screensize'));
    %% Ask to combine figure
    askcombine=input('Combine two years data? Skip (0) | Yes, but save first (1) | Yes, skip saving (2) | Skip everything and proceed next (4) ');
    %% Saving figures
    if askcombine==0 || askcombine==1
        for i=1:30
            today=char(datetime('today','Format','dd-MM-yyyy'));
            today=strcat(today,'\Adib_13_Automated_MajorEQ\');
            path1='D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\';
            mkdir(fullfile(path1,today));
            
            path2=strcat('D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Cubaan 2\Figures\',today);
            
            figname1=stnyr_sel{U};
            
            input('Hold on... Save? Enter any number ');
            
            for i=1:100
                if exist(strcat(path2,figname1,sprintf('-%d',i),'.tif'))~=2
                    figname1ext=strcat(figname1,char(sprintf('-%d',i)));
                    saveas(f1,fullfile(path2,figname1ext),'tiff');
                    disp('Figure saved')
                    break;
                end
            end
            
            askanother=input('Done? Proceed next (0) | Save another figure (1) ');
            if askanother==0 break; end
        end
    end
    %% Combine year
    if askcombine==1 || askcombine==2
        twinyear=input('Combine with what year? Previous (1) | Following (2) ');
    else
        twinyear=3;
    end
    
    if twinyear==1
        stnyr_sel{end+1}=NaN;
        for k=length(stnyr_sel)-1:-1:U+1
            stnyr_sel{k+1}=stnyr_sel{k};
        end
        stnyr_sel{U+1}=sprintf('%s%d%d',stn,yr_sel1-1,yr_sel1);
    elseif twinyear==2
        stnyr_sel{end+1}=NaN;
        for k=length(stnyr_sel)-1:-1:U+1
            stnyr_sel{k+1}=stnyr_sel{k};
        end
        stnyr_sel{U+1}=sprintf('%s%d%d',stn,yr_sel1,yr_sel1+1);
    end
    %% Clear variables
    clearvars -except U stnyr_sel
    clc
    close all
end

disp('Analysis done')