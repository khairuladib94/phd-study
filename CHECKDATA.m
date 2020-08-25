%% %Purpose
%Inspect data automatically. Search for all .mat files in the specified
%folder and plot raw and diff Z and G, nearby earthquake events as well as
%data availability percentage and total significant events (M>6.5 or dis<30km)
%then save figures in .tiff.
%% %Universal parameters/files
stn_list={'TGG';'MUT';'LGZ';'CEB';'CDO';'DAV';'GSI';'SCN';'LWA';'PTN';'MND';
    'BIK';'JYP';'PRP';'KPG';'KTB';'DAW';'LKW';'SBH';'PER';'BTG';'KTN';'TIK';
    'CHD';'CST';'ZYK';'ZGN';'MGD';'YAK';'PTK';'ASB';'TNO';'ONW';'OIS';'KUJ';
    'AMA';'HLN';'EWA';'YAP';'BCL';'HVD';'TIR';'CMB';'CKT';'TWV';'ROC';'LMT';
    'CGR';'CMD';'CAN';'MLB';'HOB';'MCQ';'DVS';'WAD';'GLY';'JRS';'TPT';'TMA';
    'ANC';'HUA';'ICA';'EUS';'SMA';'LAQ';'FYM';'ASW';'KRT';'AAB';'ILR';'ABU';
    'LAG';'ABJ';'NAB';'DES';'LSK';'MPT';'DRB';'HER'};
stn_list=string(stn_list);
yr_list=string(2005:1:2018);
filepath1='F:\Study\MAGDAS DATA\';
load VARIABLES_WORLD
mag_min=5.5;
dis_max=300;
%% %
for i=1:length(stn_list)
    for j=1:length(yr_list)
        filepath=strcat(filepath1,stn_list(i),'\');
        filename=sprintf('%3s%4sS.mat',stn_list(i),yr_list(j));
        filepathname=strcat(filepath,filename);
        filecheck=exist(filepathname,'file');
        figname1=char(sprintf('%3s%4sS_DATACHECK.tif',stn_list(i),yr_list(j)));
        figpathname=strcat(filepath,figname1);
        figcheck=exist(figpathname,'file');
        if filecheck==2 && figcheck~=2
           
            load(filename)
            
            %Remove noise/error
            G=sqrt(H.^2+D.^2);
            dnG=medfilt1(G,'omitnan');
            dnZ=medfilt1(Z,'omitnan');
            
            for k=1:5
                sigG=std(dnG,'omitnan');
                muG=mean(dnG,'omitnan');
                sigZ=std(dnZ,'omitnan');
                muZ=mean(dnZ,'omitnan');
                
                for l=1:length(dnG)
                    if dnG(l)>muG+5*sigG||dnG(l)<muG-5*sigG
                        dnG(l)=NaN;
                    end
                    if dnZ(l)>muZ+5*sigZ||dnZ(l)<muZ-5*sigZ
                        dnZ(l)=NaN;
                    end
                end
            end
            
            dataG=sum(~isnan(dnG))/length(G)*100;
            dataZ=sum(~isnan(dnZ))/length(Z)*100;
            data=sum(dataG+dataZ)/2;
            
            EQ_table=NaN(1,11);
            s=1;
            for t=1:length(EQ_WORLD)
                if EQ_WORLD(t,1)==str2num(yr_list(j))
                    EQ_table(s,:)=EQ_WORLD(t,:);
                    s=s+1;
                elseif EQ_WORLD(t,1)>str2num(yr_list(j))
                    break
                end
            end
            
            stn_latlon=[stn_MAGDAS(i,2:3)];
            t=1;
            EQ_sel=NaN(1,8);
            for s=1:size(EQ_table,1)
                if  EQ_table(s,10)>=mag_min
                    EQ_latlon=[EQ_table(s,7:8)];
                    EQ_dis=deg2km(distance('gc',stn_latlon,EQ_latlon));
                    EQ_time=datenum(EQ_table(s,1:6));
                    EQ_mag=EQ_table(s,10);
                    EQ_Ks=(10^(0.75*EQ_mag))/(EQ_dis+100);
                    EQ_depth=EQ_table(s,9);
                    if EQ_mag>=6.5 || EQ_dis<=30
                        EQ_sig=1;
                    else
                        EQ_sig=0;
                    end
                    
                    if EQ_dis<=dis_max
                        EQ_sel(t,1)=EQ_time;
                        EQ_sel(t,2)=EQ_dis;
                        EQ_sel(t,3)=EQ_mag;
                        EQ_sel(t,4)=EQ_Ks;
                        EQ_sel(t,5)=EQ_depth;
                        EQ_sel(t,6:7)=EQ_latlon;
                        EQ_sel(t,8)=EQ_sig;
                        t=t+1;
                    end
                end
            end
            
            num_sigEQ=sum(EQ_sel(:,8));
            titlestr=sprintf('Station: %s | Year: %s | Data availability: %.2f%% | Significant events: %d',stn_list(i),yr_list(j),data,num_sigEQ);
            date_vec=datetime(datevec(min(UT1m):1/86400:max(UT1m)),'Format','dd/MM/yyyy');
             
            fig=figure;
            suptitle(titlestr);
            
            subplot(4,2,[1 2])
            ZG=dnZ./dnG;
            plot(date_vec,ZG);
            ylabel('Z/G');
            xlim([date_vec(1) date_vec(end)]);
            yyaxis right
            set(gca, 'YColor', 'b')
            for u=1:size(EQ_sel,1)
                v=datetime(datevec(EQ_sel(u,1)),'Format','dd/MM/yyyy HH:mm:ss');
                if EQ_sel(u,8)==1
                    plot(v,EQ_sel(u,4),'ro','MarkerFaceColor','r','MarkerSize',5);
                else
                    plot(v,EQ_sel(u,4),'bo','MarkerFaceColor','b','MarkerSize',5);
                end
                if EQ_sel(u,4)==max(EQ_sel(:,4))
                    label=sprintf('%s M%0.1f \t%.0fKM %.0fKM AWAY',string(datetime(v,'Format','dd/MM/yyyy')),EQ_sel(u,3),EQ_sel(u,5),EQ_sel(u,2));
                    text(v,EQ_sel(u,4),label,'VerticalAlignment','top','HorizontalAlignment','left')
                end
                hold on
            end
            ylabel('K_{LS}')
            
            Z_diff=diff(dnZ); 
            Z_diff(end+1,1)=NaN;
            G_diff=diff(dnG); 
            G_diff(end+1,1)=NaN;
            
            subplot(4,2,[3 4])
            plot(date_vec,Z_diff);
            ylabel('diff(Z)');
            xlim([date_vec(1) date_vec(end)]);
            ylim([-5 5])
            
            subplot(4,2,[5 6])
            plot(date_vec,G_diff);
            ylabel('diff(G)');
            xlim([date_vec(1) date_vec(end)]);
            ylim([-5 5])
            
            n_win=1800;
            n_ovrlp=0.5*n_win;
            n_fft=n_win;
            fs=1;
            
            subplot(4,2,7)
            spectrogram(dnZ,'yaxis',hamming(n_win),n_ovrlp,n_fft,1)
            caxis([0 30])
            ylim([0 60])
            colormap jet
            colorbar('off')
            colorbar('east')
            title('Spectogram of Z')
            
            subplot(4,2,8)
            spectrogram(dnG,'yaxis',hamming(n_win),n_ovrlp,n_fft,1)
            caxis([0 30])
            ylim([0 60])
            colormap jet
            colorbar('off')
            colorbar('east')
            title('Spectogram of G')
            
            set(gcf, 'Position', get(0, 'Screensize'));
        
            figname=char(sprintf('%3s%4sS_DATACHECK',stn_list(i),yr_list(j)));
            saveas(fig,fullfile(char(filepath),figname),'tiff');
            
            display(fprintf('%s%sS FIGURE SAVED',stn_list(i),yr_list(j)));
            
            clearvars -except stn_list yr_list filepath1 i j k l EQ_WORLD stn_MAGDAS mag_min dis_max
            close all
        else
            display(fprintf('%s%sS HAS NO DATA/FIGURE ALREADY EXISTS',stn_list(i),yr_list(j)));
        end
    end
end




