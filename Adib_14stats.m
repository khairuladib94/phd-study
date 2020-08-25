%% Purpose
%To test statistical relationship of generated earthquake precursor data,
%manually
%% Read excel
xlsfilename='D:\Study\Record of appearance of earthquakes geomagnetic precursor.xlsx';
xlssheetname='Adib_14g';
[readxcel1,readxcel2,readxcel]=xlsread(xlsfilename,xlssheetname);
%% Find earthquakes with/without precursors
j=1;
k=1;
for i=2:size(readxcel,1)
    
    if strcmp(readxcel{i,15},'N/A')
        nopre(j)=i;
        j=j+1;
    else
        yespre(k)=i;
        k=k+1;
    end
end
%% Cell to array
for i=1:numel(yespre)
    yesstn1(i)=string(readxcel{yespre(i),2});
    yesmag(i)=readxcel{yespre(i),6};
    yesdep(i)=readxcel{yespre(i),7};
    yesdis(i)=readxcel{yespre(i),8};
    yesKLS(i)=readxcel{yespre(i),9};
    yestheta(i)=readxcel{yespre(i),11};
    yesregion(i)=string(readxcel{yespre(i),3});
    yesUTC(i)=string(readxcel{yespre(i),5});
    yesLT(i)=string(readxcel{yespre(i),16});
    yesZG(i)=readxcel{yespre(i),17};
    yesZ(i)=readxcel{yespre(i),18};
    yesG(i)=readxcel{yespre(i),19};
    yeslead(i)=readxcel{yespre(i),20};
    yesfreq(i)=readxcel{yespre(i),21};
    yesorder(i)=readxcel{yespre(i),24};
    yeslat(i)=readxcel{yespre(i),12};
    yeslon(i)=readxcel{yespre(i),13};
    yesadjZG_1(i)=readxcel{yespre(i),29};
    yesadjZG_2(i)=readxcel{yespre(i),30};
end
yesstn2=char(yesstn1); yesstn=string(yesstn2(1,1:3,:));
yesUTC1=char(yesUTC); 
for i=1:size(yesUTC1,3) yesyear(i)=str2double(strcat('20',yesUTC1(1,8:9,i))); end
clearvars yesstn1 yesstn2 yesUTC1 

for i=1:numel(nopre)
    nomag(i)=readxcel{nopre(i),6};
    nodep(i)=readxcel{nopre(i),7};
    nodis(i)=readxcel{nopre(i),8};
    noKLS(i)=readxcel{nopre(i),9};
    notheta(i)=readxcel{nopre(i),11};
    noregion(i)=string(readxcel{nopre(i),3});
    noUTC(i)=string(readxcel{nopre(i),5});
    nolat(i)=readxcel{nopre(i),12};
    nolon(i)=readxcel{nopre(i),13};
end
%% #1: Pie chart of occurences of earthquake precursors
PIEoccall=figure(1);
label={sprintf('With precursor (%d)',numel(yespre)),sprintf('Without precursor (%d)',numel(nopre))};
pie1=pie([numel(yespre),numel(nopre)],label);
pie1(1).FaceColor='b'; pie1(3).FaceColor='r';
title('Pie chart of occurences of earthquake precursors');

% set(gcf, 'Position', get(0, 'Screensize'));
%% #2: Bar chart of earthquake precursor occurences versus magnitude
k=1;
for i=5.5:0.5:9.0
    yespremag(k,1)=i+0.25;
    yespremag(k,2)=sum(yesmag>=i & yesmag<i+0.5);
    nopremag(k,1)=i+0.25;
    nopremag(k,2)=sum(nomag>=i & nomag<i+0.5);
    k=k+1;
end

BARoccVSmag=figure(2);
stitle=suptitle('Bar chart of earthquake precursor occurences versus magnitude');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);


subplot(2,2,[3 4])
k=1;
for i=5.5:0.1:9.0
    yespremag1(k,1)=i;
    yespremag1(k,2)=sum(yesmag==i);
    nopremag1(k,1)=i;
    nopremag1(k,2)=sum(nomag==i);
    k=k+1;
end

plot(yespremag1(:,1),yespremag1(:,2),'bo-');
hold on
plot(nopremag1(:,1),nopremag1(:,2),'ro-');
hold off
set(gca,'XLim',[5.5 9.0],'XTick',5.5:0.2:9.0);
posi=get(gca,'Position');
set(gca,'Position',[posi(1) 0.07 posi(3) 0.2])
xlabel('Magnitude');
ylabel('Occurences')

subplot(2,2,1)
bar(yespremag(:,1),yespremag(:,2),'b');
ylim([0 15]); xlim([5.5 9.5]);
title('Earthquakes with precursor');
xlabel('Magnitude')
ylabel('Occurences')
posi=get(gca,'Position');
set(gca,'Position',[posi(1) posi(2)-0.23 posi(3) posi(4)+0.23]);

subplot(2,2,2)
bar(nopremag(:,1),nopremag(:,2),'r');
ylim([0 max(max(yespremag(:,2),nopremag(:,2)))+1]); xlim([5.5 9.5]);
title('Earthquakes without precursor');
xlabel('Magnitude')
ylabel('Occurences')
posi=get(gca,'Position');
set(gca,'Position',[posi(1) posi(2)-0.23 posi(3) posi(4)+0.23]);

% set(gcf, 'Position', get(0, 'Screensize'));
%% #3: Bar chart of earthquake precursor occurences versus depth

k=1;
for i=0:20:180
    yespredep(k,1)=i+10;
    yespredep(k,2)=sum(yesdep>=i & yesdep<i+20);
    nopredep(k,1)=i+10;
    nopredep(k,2)=sum(nodep>=i & nodep<i+20);
    k=k+1;
end

BARoccVSdep=figure(3)
stitle=suptitle('Bar chart of earthquake precursor occurences versus depth');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

subplot(2,2,[3 4])
plot(yesdep,zeros(1,numel(yesdep)),'bo','MarkerSize',8);
hold on
plot(nodep,zeros(1,numel(nodep)),'ro','MarkerSize',8);
hold off
set(gca,'DataAspectRatio',[1 1 1],'XLim',[0 200],...
    'YLim',[0 eps],'Color','none','XTick',0:20:200);
posi=get(gca,'Position');
set(gca,'Position',[posi(1) 0.07 posi(3) 0.1])
xlabel('Depth (km)');

subplot(2,2,1)
bar(yespredep(:,1),yespredep(:,2),'b');
ylim([0 max(max(yespredep(:,2),nopredep(:,2)))+1]); xlim([0 200]);
title('Earthquakes with precursor');
xlabel('Depth (km)')
ylabel('Occurences')
xticks(0:20:200);
posi=get(gca,'Position');
set(gca,'Position',[posi(1) posi(2)-0.35 posi(3) posi(4)+0.35])


subplot(2,2,2)
bar(nopredep(:,1),nopredep(:,2),'r');
ylim([0 max(max(yespredep(:,2),nopredep(:,2)))+1]); xlim([0 200]);
title('Earthquakes without precursor');
xlabel('Depth (km)')
ylabel('Occurences')
xticks(0:20:200);
posi=get(gca,'Position');
set(gca,'Position',[posi(1) posi(2)-0.35 posi(3) posi(4)+0.35])

% set(gcf, 'Position', get(0, 'Screensize'));
%% #4: Bar chart of earthquake precursor occurences versus epicentral distance
k=1;
for i=0:20:200
    yespredis(k,1)=i+10;
    yespredis(k,2)=sum(yesdis>=i & yesdis<i+20);
    nopredis(k,1)=i+10;
    nopredis(k,2)=sum(nodis>=i & nodis<i+20);
    k=k+1;
end

BARoccVSdis=figure(4)
stitle=suptitle('Bar chart of earthquake precursor occurences versus epicentral distance');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

subplot(2,2,[3 4])
plot(yesdis,zeros(1,numel(yesdis)),'bo','MarkerSize',8);
hold on
plot(nodis,zeros(1,numel(nodis)),'ro','MarkerSize',8);
hold off
set(gca,'DataAspectRatio',[1 1 1],'XLim',[0 200],...
    'YLim',[0 eps],'Color','none','XTick',0:20:200);
posi=get(gca,'Position');
set(gca,'Position',[posi(1) 0.07 posi(3) 0.1])
xlabel('Distance (km)');

subplot(2,2,1)
bar(yespredis(:,1),yespredis(:,2),'b');
ylim([0 max(max(yespredis(:,2),nopredis(:,2)))+1]); xlim([0 200]);
title('Earthquakes with precursor');
xlabel('Distance (km)')
ylabel('Occurences')
xticks(0:20:200);
posi=get(gca,'Position');
set(gca,'Position',[posi(1) posi(2)-0.35 posi(3) posi(4)+0.35])


subplot(2,2,2)
bar(nopredis(:,1),nopredis(:,2),'r');
ylim([0 max(max(yespredis(:,2),nopredis(:,2)))+1]); xlim([0 200]);
title('Earthquakes without precursor');
xlabel('Distance (km)')
ylabel('Occurences')
xticks(0:20:200);
posi=get(gca,'Position');
set(gca,'Position',[posi(1) posi(2)-0.35 posi(3) posi(4)+0.35])

% set(gcf, 'Position', get(0, 'Screensize'));
%% #5: Bar chart of earthquake precursor occurences versus seismicity index, KLS
k=1;
for i=0:500:max([max(yesKLS) max(noKLS)])
    yespreKLS(k,1)=i+250;
    yespreKLS(k,2)=sum(yesKLS>=i & yesKLS<i+500);
    nopreKLS(k,1)=i+250;
    nopreKLS(k,2)=sum(noKLS>=i & noKLS<i+500);
    k=k+1;
end

BARoccVSkls=figure(5)
stitle=suptitle('Bar chart of earthquake precursor occurences versus seismicity index, K_{LS}');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);


subplot(2,2,[3 4])
plot(yesKLS,zeros(1,numel(yesKLS)),'bo','MarkerSize',8);
hold on
plot(noKLS,zeros(1,numel(noKLS)),'ro','MarkerSize',8);
hold off
set(gca,'DataAspectRatio',[1 1 1],'XLim',[0 max([max(yesKLS) max(noKLS)])],...
    'YLim',[0 eps],'Color','none');
posi=get(gca,'Position');
set(gca,'Position',[posi(1) 0.07 posi(3) 0.1])
xlabel('K_{LS}')


subplot(2,2,1)
bar(yespreKLS(:,1),yespreKLS(:,2),'b');
ylim([0 max(max(yespreKLS(:,2),nopreKLS(:,2)))+1]); xlim([0 max([max(yesKLS) max(noKLS)])]);
title('Earthquakes with precursor');
xlabel('K_{LS}')
ylabel('Occurences')
posi=get(gca,'Position');
set(gca,'Position',[posi(1) posi(2)-0.35 posi(3) posi(4)+0.35])

subplot(2,2,2)
bar(nopreKLS(:,1),nopreKLS(:,2),'r');
ylim([0 max(max(yespreKLS(:,2),nopreKLS(:,2)))+1]); xlim([0 max([max(yesKLS) max(noKLS)])]);
title('Earthquakes without precursor');
xlabel('Distance (km)')
ylabel('Occurences')
posi=get(gca,'Position');
set(gca,'Position',[posi(1) posi(2)-0.35 posi(3) posi(4)+0.35])

% set(gcf, 'Position', get(0, 'Screensize'));
%% #6: Bar chart of earthquake precursor occurences versus azimuthal angle
k=1;
for i=0:45:315
    yespretheta(k,1)=i+22.5;
    yespretheta(k,2)=sum(yestheta>=i & yestheta<i+45);
    nopretheta(k,1)=i+22.5;
    nopretheta(k,2)=sum(notheta>=i & notheta<i+45);
    k=k+1;
end

BARoccVStheta=figure(6)
stitle=suptitle('Bar chart of earthquake precursor occurences versus azimuthal angle');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);


subplot(2,2,[3 4])
plot(yestheta,zeros(1,numel(yestheta)),'bo','MarkerSize',8);
hold on
plot(notheta,zeros(1,numel(notheta)),'ro','MarkerSize',8);
hold off
set(gca,'DataAspectRatio',[1 1 1],'XLim',[0 360],...
    'YLim',[0 eps],'Color','none','XTick',0:45:360);
posi=get(gca,'Position');
set(gca,'Position',[posi(1) 0.07 posi(3) 0.1])
xlabel('Azimuthal angle (\circ)');


subplot(2,2,1)
bar(yespretheta(:,1),yespretheta(:,2),'b');
ylim([0 max(max(yespretheta(:,2),nopretheta(:,2)))+1]); xlim([0 360]);
title('Earthquakes with precursor');
xlabel('Azimuthal angle (\circ)');
ylabel('Occurences');
xticks(0:45:360)
posi=get(gca,'Position');
set(gca,'Position',[posi(1) posi(2)-0.35 posi(3) posi(4)+0.35])

subplot(2,2,2)
bar(nopretheta(:,1),nopretheta(:,2),'r');
ylim([0 max(max(yespretheta(:,2),nopretheta(:,2)))+1]); xlim([0 360]);
title('Earthquakes without precursor');
xlabel('Azimuthal angle (\circ)')
ylabel('Occurences');
xticks(0:45:360);
posi=get(gca,'Position');
set(gca,'Position',[posi(1) posi(2)-0.35 posi(3) posi(4)+0.35])

% set(gcf, 'Position', get(0, 'Screensize'));
%% #7: Bar chart of earthquake precursor occurences based on region and year

k=1;
regionarray=unique(horzcat(yesregion,noregion));
for i=1:numel(regionarray)
    yespreregion(k,1)=regionarray(i);
    yespreregion(k,2)=sum(yesregion==regionarray(i));
    nopreregion(k,1)=regionarray(i);
    nopreregion(k,2)=sum(noregion==regionarray(i));
    k=k+1;
end
k=1;
charyesUTC=char(yesUTC); charnoUTC=char(noUTC);
for i=2005:2018
    j=char(string(i));
    yespreyear(k,1)=i;
    yespreyear(k,2)=sum(string(charyesUTC(:,8:9,:))==string(j(3:4)));
    nopreyear(k,1)=i;
    nopreyear(k,2)=sum(string(charnoUTC(:,8:9,:))==string(j(3:4)));
    k=k+1;
end

BARoccVSregionyear=figure(7)
stitle=suptitle('Bar chart of earthquake precursor occurences based on region and year');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

subplot(2,1,1)
bar7a=bar(categorical(yespreregion(:,1)),horzcat(str2double(yespreregion(:,2)),str2double(nopreregion(:,2))),'stacked');
bar7a(1).FaceColor='blue'; bar7a(2).FaceColor='red';
xlabel('Region');
ylabel('Occurences');
legend({'With precursor','Without precursor'})

yearmin=yespreyear(min([min(find(yespreyear(:,2)~=0)) min(find(nopreyear(:,2)~=0))]),1);
yearmax=yespreyear(max([max(find(yespreyear(:,2)~=0)) max(find(nopreyear(:,2)~=0))]),1);

subplot(2,1,2)
bar7b=bar(yespreyear(:,1),horzcat(yespreyear(:,2),nopreyear(:,2)),'stacked');
bar7b(1).FaceColor='blue'; bar7b(2).FaceColor='red';
xlabel('Year');
ylabel('Occurences');
legend({'With precursor','Without precursor'})
xlim([yearmin-.5 yearmax+.5])

% set(gcf, 'Position', get(0, 'Screensize'));
%% #8: Scatter plot of amplitude of Z/G versus magnitude, depth, distance and KLS
SCATzgVSmagdepdisKLS=figure(8);

stitle=suptitle('Scatter plot of amplitude of Z/G versus magnitude, depth, distance and K_{LS}');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

subplot(2,2,1)
[FITmagZG,GOFmagZG]=fit(yesmag',yesZG','poly1');
magZG_fit=FITmagZG.p1*yesmag+FITmagZG.p2;
scatter(yesmag,yesZG,'bo','HandleVisibility','off');
hold on
plot(yesmag,magZG_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFmagZG.rsquare,GOFmagZG.rmse));
text(max(yesmag),max(yesZG),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Z/G');
xlabel('Magnitude')

subplot(2,2,2)
[FITdepZG,GOFdepZG]=fit(yesdep',yesZG','poly1');
depZG_fit=FITdepZG.p1*yesdep+FITdepZG.p2;
scatter(yesdep,yesZG,'bo','HandleVisibility','off');
hold on
plot(yesdep,depZG_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFdepZG.rsquare,GOFdepZG.rmse));
text(max(yesdep),max(yesZG),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Z/G');
xlabel('Depth (km)')

subplot(2,2,3)
[FITdisZG,GOFdisZG]=fit(yesdis',yesZG','poly1');
disZG_fit=FITdisZG.p1*yesdis+FITdisZG.p2;
scatter(yesdis,yesZG,'bo','HandleVisibility','off');
hold on
plot(yesdis,disZG_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFdisZG.rsquare,GOFdisZG.rmse));
text(max(yesdis),max(yesZG),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Z/G');
xlabel('Distance (km)')

subplot(2,2,4)
[FITKLSZG,GOFKLSZG]=fit(yesKLS',yesZG','poly1');
KLSZG_fit=FITKLSZG.p1*yesKLS+FITKLSZG.p2;
scatter(yesKLS,yesZG,'bo','HandleVisibility','off');
hold on
plot(yesKLS,KLSZG_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFKLSZG.rsquare,GOFKLSZG.rmse));
text(max(yesKLS),max(yesZG),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Z/G');
xlabel('K_{LS}')

% set(gcf, 'Position', get(0, 'Screensize'));
%% #9: Scatter plot of amplitude of Z/G versus latitude and longitude of earthquake
SCATzgVSlatlon=figure(89);

stitle=suptitle('Scatter plot of amplitude of Z/G versus latitude and longitude of earthquake');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

nameregion=unique(yesregion);
colpal=jet(numel(nameregion));

subplot(2,2,1)
for i=1:numel(nameregion)
    hold on
    plot(yeslat(find(strcmp(yesregion,nameregion(i)))),yesZG(find(strcmp(yesregion,nameregion(i)))),'o','Color',colpal(i,:));
    hold off
end
set(gca,'YAxisLocation','origin');
xlim([-90 90])
ylabel('Z/G');
xlabel('Latitude')

subplot(2,2,2)
for i=1:numel(nameregion)
    hold on
    plot(yeslon(find(strcmp(yesregion,nameregion(i)))),yesZG(find(strcmp(yesregion,nameregion(i)))),'o','Color',colpal(i,:));
    hold off
end
set(gca,'YAxisLocation','origin');
xlim([-180 180])
legend(nameregion,'Location','northwest')
ylabel('Z/G');
xlabel('Longitude')

subplot(2,2,[3 4]);
for i=1:numel(nameregion)
    hold on
    scatter3(yeslon(find(strcmp(yesregion,nameregion(i)))),yeslat(find(strcmp(yesregion,nameregion(i)))),yesZG(find(strcmp(yesregion,nameregion(i)))),'o','MarkerEdgeColor',colpal(i,:));
    hold off
end
ylim([-90 90]); xlim([-180 180]);
xlabel('Longitude'); ylabel('Latitude'); zlabel('Z/G');
grid('on');
view(5,35)

% set(gcf, 'Position', get(0, 'Screensize'));
%% #10 & #11: Scatter plot of amplitude of Z/G versus regional magnitude

nameregion=unique(yesregion);
colpal=jet(numel(nameregion));

SCATzgVSregionalmag1=figure(10)

stitle=suptitle('Scatter plot of amplitude of Z/G versus regional magnitude');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=1:4
    subplot(2,2,i)
    yesmag_region=yesmag(find(strcmp(yesregion,nameregion(i))));
    yesZG_region=yesZG(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesmag_region,yesZG_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesmag_region)>1
        [FIT,GOF]=fit(yesmag_region',yesZG_region','poly1');
        [Rval,Pval]=corrcoef(yesmag_region,yesZG_region);
        fitline=FIT.p1*yesmag_region+FIT.p2;
        plot(yesmag_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesmag),max(yesZG),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([min(yesZG)-1 max(yesZG)+1])
    xlim([min(yesmag)-0.5 max(yesmag)+0.5])
    ylabel('Z/G');
    xlabel('Magnitude')
    title(sprintf('%s',nameregion(i)))
end
hold off
ylim([min(yesZG)-1 max(yesZG)+1])
xlim([min(yesmag)-0.5 max(yesmag)+0.5])
ylabel('Z/G');
xlabel('Magnitude')


SCATzgVSregionalmag2=figure(11)

stitle=suptitle('Scatter plot of amplitude of Z/G versus regional magnitude');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=5:8
    subplot(2,2,i-4)
    yesmag_region=yesmag(find(strcmp(yesregion,nameregion(i))));
    yesZG_region=yesZG(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesmag_region,yesZG_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesmag_region)>1
        [FIT,GOF]=fit(yesmag_region',yesZG_region','poly1');
        [Rval,Pval]=corrcoef(yesmag_region,yesZG_region);
        fitline=FIT.p1*yesmag_region+FIT.p2;
        plot(yesmag_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesmag),max(yesZG),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([min(yesZG)-1 max(yesZG)+1])
    xlim([min(yesmag)-0.5 max(yesmag)+0.5])
    ylabel('Z/G');
    xlabel('Magnitude')
    title(sprintf('%s',nameregion(i)))
end

% set(gcf, 'Position', get(0, 'Screensize'));
%% #12 & #13: Scatter plot of amplitude of Z/G versus regional distance

nameregion=unique(yesregion);
colpal=jet(numel(nameregion));

SCATzgVSregionaldis1=figure(12)

stitle=suptitle('Scatter plot of amplitude of Z/G versus regional distance');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=1:4
    subplot(2,2,i)
    yesdis_region=yesdis(find(strcmp(yesregion,nameregion(i))));
    yesZG_region=yesZG(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesdis_region,yesZG_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesdis_region)>1
        [FIT,GOF]=fit(yesdis_region',yesZG_region','poly1');
        [Rval,Pval]=corrcoef(yesdis_region,yesZG_region);
        fitline=FIT.p1*yesdis_region+FIT.p2;
        plot(yesdis_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesdis),max(yesZG),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([min(yesZG)-1 max(yesZG)+1])
    xlim([min(yesdis)-10 max(yesdis)+10])
    ylabel('Z/G');
    xlabel('Distance (km)')
    title(sprintf('%s',nameregion(i)))
end

SCATzgVSregionaldis2=figure(13)

stitle=suptitle('Scatter plot of amplitude of Z/G versus regional distance');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=5:8
    subplot(2,2,i-4)
    yesdis_region=yesdis(find(strcmp(yesregion,nameregion(i))));
    yesZG_region=yesZG(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesdis_region,yesZG_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesdis_region)>1
        [FIT,GOF]=fit(yesdis_region',yesZG_region','poly1');
        [Rval,Pval]=corrcoef(yesdis_region,yesZG_region);
        fitline=FIT.p1*yesdis_region+FIT.p2;
        plot(yesdis_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesdis),max(yesZG),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([min(yesZG)-1 max(yesZG)+1])
    xlim([min(yesdis)-10 max(yesdis)+10])
    ylabel('Z/G');
    xlabel('Distance (km)')
    title(sprintf('%s',nameregion(i)))
end

% set(gcf, 'Position', get(0, 'Screensize'));
%% #14 & #15: Scatter plot of amplitude of Z/G versus regional depth

nameregion=unique(yesregion);
colpal=jet(numel(nameregion));

SCATzgVSregionaldep1=figure(14)

stitle=suptitle('Scatter plot of amplitude of Z/G versus regional depth');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=1:4
    subplot(2,2,i)
    yesdep_region=yesdep(find(strcmp(yesregion,nameregion(i))));
    yesZG_region=yesZG(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesdep_region,yesZG_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesdep_region)>1
        [FIT,GOF]=fit(yesdep_region',yesZG_region','poly1');
        [Rval,Pval]=corrcoef(yesdep_region,yesZG_region);
        fitline=FIT.p1*yesdep_region+FIT.p2;
        plot(yesdep_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesdep),max(yesZG),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([min(yesZG)-1 max(yesZG)+1])
    xlim([min(yesdep)-10 max(yesdep)+10])
    ylabel('Z/G');
    xlabel('Depth (km)')
    title(sprintf('%s',nameregion(i)))
end

SCATzgVSregionaldep2=figure(15)

stitle=suptitle('Scatter plot of amplitude of Z/G versus regional depth');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=5:8
    subplot(2,2,i-4)
    yesdep_region=yesdep(find(strcmp(yesregion,nameregion(i))));
    yesZG_region=yesZG(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesdep_region,yesZG_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesdep_region)>1
        [FIT,GOF]=fit(yesdep_region',yesZG_region','poly1');
        [Rval,Pval]=corrcoef(yesdep_region,yesZG_region);
        fitline=FIT.p1*yesdep_region+FIT.p2;
        plot(yesdep_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesdep),max(yesZG),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([min(yesZG)-1 max(yesZG)+1])
    xlim([min(yesdep)-10 max(yesdep)+10])
    ylabel('Z/G');
    xlabel('Depth (km)')
    title(sprintf('%s',nameregion(i)))
end

% set(gcf, 'Position', get(0, 'Screensize'));
%% #16: Scatter plot of effective frequency versus magnitude, depth, distance and KLS
SCATfreqVSmagdepdisKLS=figure(16);

stitle=suptitle('Scatter plot of effective frequency versus magnitude, depth, distance and K_{LS}');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

subplot(2,2,1)
[FITmagfreq,GOFmagfreq]=fit(yesmag',yesfreq','poly1');
magfreq_fit=FITmagfreq.p1*yesmag+FITmagfreq.p2;
scatter(yesmag,yesfreq,'bo','HandleVisibility','off');
hold on
plot(yesmag,magfreq_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFmagfreq.rsquare,GOFmagfreq.rmse));
text(max(yesmag),max(yesfreq),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Frequency (Hz)');
xlabel('Magnitude')

subplot(2,2,2)
[FITdepfreq,GOFdepfreq]=fit(yesdep',yesfreq','poly1');
depfreq_fit=FITdepfreq.p1*yesdep+FITdepfreq.p2;
scatter(yesdep,yesfreq,'bo','HandleVisibility','off');
hold on
plot(yesdep,depfreq_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFdepfreq.rsquare,GOFdepfreq.rmse));
text(max(yesdep),max(yesfreq),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Frequency (Hz)');
xlabel('Depth (km)')

subplot(2,2,3)
[FITdisfreq,GOFdisfreq]=fit(yesdis',yesfreq','poly1');
disfreq_fit=FITdisfreq.p1*yesdis+FITdisfreq.p2;
scatter(yesdis,yesfreq,'bo','HandleVisibility','off');
hold on
plot(yesdis,disfreq_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFdisfreq.rsquare,GOFdisfreq.rmse));
text(max(yesdis),max(yesfreq),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Frequency (Hz)');
xlabel('Distance (km)')

subplot(2,2,4)
[FITKLSfreq,GOFKLSfreq]=fit(yesKLS',yesfreq','poly1');
KLSfreq_fit=FITKLSfreq.p1*yesKLS+FITKLSfreq.p2;
scatter(yesKLS,yesfreq,'bo','HandleVisibility','off');
hold on
plot(yesKLS,KLSfreq_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFKLSfreq.rsquare,GOFKLSfreq.rmse));
text(max(yesKLS),max(yesfreq),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Frequency (Hz)');
xlabel('K_{LS}')

% set(gcf, 'Position', get(0, 'Screensize'));
%% #17: Scatter plot of effective frequency versus latitude and longitude of earthquake
SCATfreqVSlatlon=figure(17);

stitle=suptitle('Scatter plot of effective frequency versus latitude and longitude of earthquake');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

nameregion=unique(yesregion);
colpal=jet(numel(nameregion));

subplot(2,2,1)
for i=1:numel(nameregion)
    hold on
    plot(yeslat(find(strcmp(yesregion,nameregion(i)))),yesfreq(find(strcmp(yesregion,nameregion(i)))),'o','Color',colpal(i,:));
    hold off
end
set(gca,'YAxisLocation','origin');
xlim([-90 90])
ylabel('Frequency (Hz)');
xlabel('Latitude')

subplot(2,2,2)
for i=1:numel(nameregion)
    hold on
    plot(yeslon(find(strcmp(yesregion,nameregion(i)))),yesfreq(find(strcmp(yesregion,nameregion(i)))),'o','Color',colpal(i,:));
    hold off
end
set(gca,'YAxisLocation','origin');
xlim([-180 180])
legend(nameregion,'Location','northwest')
ylabel('Frequency (Hz)');
xlabel('Longitude')

subplot(2,2,[3 4]);
for i=1:numel(nameregion)
    hold on
    scatter3(yeslon(find(strcmp(yesregion,nameregion(i)))),yeslat(find(strcmp(yesregion,nameregion(i)))),yesfreq(find(strcmp(yesregion,nameregion(i)))),'o','MarkerEdgeColor',colpal(i,:));
    hold off
end
ylim([-90 90]); xlim([-180 180]);
xlabel('Longitude'); ylabel('Latitude'); zlabel('Frequency (Hz)');
grid('on');
view(5,35)

% set(gcf, 'Position', get(0, 'Screensize'));
%% #18 & #19: Scatter plot of effective frequency versus regional magnitude

nameregion=unique(yesregion);
colpal=jet(numel(nameregion));

SCATfreqVSregionalmag1=figure(18)

stitle=suptitle('Scatter plot of effective frequency versus regional magnitude');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=1:4
    subplot(2,2,i)
    yesmag_region=yesmag(find(strcmp(yesregion,nameregion(i))));
    yesfreq_region=yesfreq(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesmag_region,yesfreq_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesmag_region)>1
        [FIT,GOF]=fit(yesmag_region',yesfreq_region','poly1');
        [Rval,Pval]=corrcoef(yesmag_region,yesfreq_region);
        fitline=FIT.p1*yesmag_region+FIT.p2;
        plot(yesmag_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesmag),max(yesfreq),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([0.01 0.1])
    xlim([min(yesmag)-0.5 max(yesmag)+0.5])
    ylabel('Frequency (Hz)');
    xlabel('Magnitude')
    title(sprintf('%s',nameregion(i)))
end

SCATfreqVSregionalmag2=figure(19)

stitle=suptitle('Scatter plot of effective frequency versus regional magnitude');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=5:8
    subplot(2,2,i-4)
    yesmag_region=yesmag(find(strcmp(yesregion,nameregion(i))));
    yesfreq_region=yesfreq(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesmag_region,yesfreq_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesmag_region)>1
        [FIT,GOF]=fit(yesmag_region',yesfreq_region','poly1');
        [Rval,Pval]=corrcoef(yesmag_region,yesfreq_region);
        fitline=FIT.p1*yesmag_region+FIT.p2;
        plot(yesmag_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesmag),max(yesfreq),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([0.01 0.1])
    xlim([min(yesmag)-0.5 max(yesmag)+0.5])
    ylabel('Frequency (Hz)');
    xlabel('Magnitude')
    title(sprintf('%s',nameregion(i)))
end

% set(gcf, 'Position', get(0, 'Screensize'));
%% #20 & #21: Scatter plot of effective frequency versus regional depth

nameregion=unique(yesregion);
colpal=jet(numel(nameregion));

SCATfreqVSregionaldep1=figure(20)

stitle=suptitle('Scatter plot of effective frequency versus regional depth');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=1:4
    subplot(2,2,i)
    yesdep_region=yesdep(find(strcmp(yesregion,nameregion(i))));
    yesfreq_region=yesfreq(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesdep_region,yesfreq_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesdep_region)>1
        [FIT,GOF]=fit(yesdep_region',yesfreq_region','poly1');
        [Rval,Pval]=corrcoef(yesdep_region,yesfreq_region);
        fitline=FIT.p1*yesdep_region+FIT.p2;
        plot(yesdep_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesdep),max(yesfreq),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([0.01 0.1])
    xlim([min(yesdep)-10 max(yesdep)+10])
    ylabel('Frequency (Hz)');
    xlabel('Depth (km)')
    title(sprintf('%s',nameregion(i)))
end

SCATfreqVSregionaldep2=figure(21)

stitle=suptitle('Scatter plot of effective frequency versus regional depth');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=5:8
    subplot(2,2,i-4)
    yesdep_region=yesdep(find(strcmp(yesregion,nameregion(i))));
    yesfreq_region=yesfreq(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesdep_region,yesfreq_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesdep_region)>1
        [FIT,GOF]=fit(yesdep_region',yesfreq_region','poly1');
        [Rval,Pval]=corrcoef(yesdep_region,yesfreq_region);
        fitline=FIT.p1*yesdep_region+FIT.p2;
        plot(yesdep_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesdep),max(yesfreq),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([0.01 0.1])
    xlim([min(yesdep)-10 max(yesdep)+10])
    ylabel('Frequency (Hz)');
    xlabel('Depth (km)')
    title(sprintf('%s',nameregion(i)))
end

% set(gcf, 'Position', get(0, 'Screensize'));
%% #22 & #23: Scatter plot of effective frequency versus regional distance

nameregion=unique(yesregion);
colpal=jet(numel(nameregion));

SCATfreqVSregionaldis1=figure(22)

stitle=suptitle('Scatter plot of effective frequency versus regional distance');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=1:4
    subplot(2,2,i)
    yesdis_region=yesdis(find(strcmp(yesregion,nameregion(i))));
    yesfreq_region=yesfreq(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesdis_region,yesfreq_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesdis_region)>1
        [FIT,GOF]=fit(yesdis_region',yesfreq_region','poly1');
        [Rval,Pval]=corrcoef(yesdis_region,yesfreq_region);
        fitline=FIT.p1*yesdis_region+FIT.p2;
        plot(yesdis_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesdis),max(yesfreq),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([0.01 0.1])
    xlim([min(yesdis)-10 max(yesdis)+10])
    ylabel('Frequency (Hz)');
    xlabel('Distance (km)')
    title(sprintf('%s',nameregion(i)))
end

SCATfreqVSregionaldis2=figure(23)

stitle=suptitle('Scatter plot of effective frequency versus regional distance');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=5:8
    subplot(2,2,i-4)
    yesdis_region=yesdis(find(strcmp(yesregion,nameregion(i))));
    yesfreq_region=yesfreq(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesdis_region,yesfreq_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesdis_region)>1
        [FIT,GOF]=fit(yesdis_region',yesfreq_region','poly1');
        [Rval,Pval]=corrcoef(yesdis_region,yesfreq_region);
        fitline=FIT.p1*yesdis_region+FIT.p2;
        plot(yesdis_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesdis),max(yesfreq),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([0.01 0.1])
    xlim([min(yesdis)-10 max(yesdis)+10])
    ylabel('Frequency (Hz)');
    xlabel('Distance (km)')
    title(sprintf('%s',nameregion(i)))
end

% set(gcf, 'Position', get(0, 'Screensize'));
%% #24: Scatter plot of lead time versus magnitude, depth, distance and KLS

SCATleadVSmagdepdisKLS=figure(24);

stitle=suptitle('Scatter plot of lead time versus magnitude, depth, distance and K_{LS}');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

subplot(2,2,1)
[FITmaglead,GOFmaglead]=fit(yesmag',yeslead','poly1');
maglead_fit=FITmaglead.p1*yesmag+FITmaglead.p2;
scatter(yesmag,yeslead,'bo','HandleVisibility','off');
hold on
plot(yesmag,maglead_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFmaglead.rsquare,GOFmaglead.rmse));
text(max(yesmag),max(yeslead),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Lead time (days)');
xlabel('Magnitude')

subplot(2,2,2)
[FITdeplead,GOFdeplead]=fit(yesdep',yeslead','poly1');
deplead_fit=FITdeplead.p1*yesdep+FITdeplead.p2;
scatter(yesdep,yeslead,'bo','HandleVisibility','off');
hold on
plot(yesdep,deplead_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFdeplead.rsquare,GOFdeplead.rmse));
text(max(yesdep),max(yeslead),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Lead time (days)');
xlabel('Depth (km)')

subplot(2,2,3)
[FITdislead,GOFdislead]=fit(yesdis',yeslead','poly1');
dislead_fit=FITdislead.p1*yesdis+FITdislead.p2;
scatter(yesdis,yeslead,'bo','HandleVisibility','off');
hold on
plot(yesdis,dislead_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFdislead.rsquare,GOFdislead.rmse));
text(max(yesdis),max(yeslead),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Lead time (days)');
xlabel('Distance (km)')

subplot(2,2,4)
[FITKLSlead,GOFKLSlead]=fit(yesKLS',yeslead','poly1');
KLSlead_fit=FITKLSlead.p1*yesKLS+FITKLSlead.p2;
scatter(yesKLS,yeslead,'bo','HandleVisibility','off');
hold on
plot(yesKLS,KLSlead_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFKLSlead.rsquare,GOFKLSlead.rmse));
text(max(yesKLS),max(yeslead),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Lead time (days)');
xlabel('K_{LS}')

% set(gcf, 'Position', get(0, 'Screensize'));
%% #25: Scatter plot of lead time versus latitude and longitude of earthquake

SCATleadVSlatlon=figure(25);

stitle=suptitle('Scatter plot of lead time versus latitude and longitude of earthquake');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

nameregion=unique(yesregion);
colpal=jet(numel(nameregion));

subplot(2,2,1)
for i=1:numel(nameregion)
    hold on
    plot(yeslat(find(strcmp(yesregion,nameregion(i)))),yeslead(find(strcmp(yesregion,nameregion(i)))),'o','Color',colpal(i,:));
    hold off
end
set(gca,'YAxisLocation','origin');
xlim([-90 90])
ylabel('Lead time (days)');
xlabel('Latitude')

subplot(2,2,2)
for i=1:numel(nameregion)
    hold on
    plot(yeslon(find(strcmp(yesregion,nameregion(i)))),yeslead(find(strcmp(yesregion,nameregion(i)))),'o','Color',colpal(i,:));
    hold off
end
set(gca,'YAxisLocation','origin');
xlim([-180 180])
legend(nameregion,'Location','northwest')
ylabel('Lead time (days)');
xlabel('Longitude')

subplot(2,2,[3 4]);
for i=1:numel(nameregion)
    hold on
    scatter3(yeslon(find(strcmp(yesregion,nameregion(i)))),yeslat(find(strcmp(yesregion,nameregion(i)))),yeslead(find(strcmp(yesregion,nameregion(i)))),'o','MarkerEdgeColor',colpal(i,:));
    hold off
end
ylim([-90 90]); xlim([-180 180]);
xlabel('Longitude'); ylabel('Latitude'); zlabel('Lead time (days)');
grid('on');
view(5,35)

% set(gcf, 'Position', get(0, 'Screensize'));
%% #26 & #27: Scatter plot of lead time versus regional magnitude

nameregion=unique(yesregion);
colpal=jet(numel(nameregion));

SCATleadVSregionalmag1=figure(26)

stitle=suptitle('Scatter plot of lead time versus regional magnitude');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=1:4
    subplot(2,2,i)
    yesmag_region=yesmag(find(strcmp(yesregion,nameregion(i))));
    yeslead_region=yeslead(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesmag_region,yeslead_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesmag_region)>1
        [FIT,GOF]=fit(yesmag_region',yeslead_region','poly1');
        [Rval,Pval]=corrcoef(yesmag_region,yeslead_region);
        fitline=FIT.p1*yesmag_region+FIT.p2;
        plot(yesmag_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesmag),max(yeslead),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([min(yeslead) max(yeslead)])
    xlim([min(yesmag)-0.5 max(yesmag)+0.5])
    ylabel('Lead time (days)');
    xlabel('Magnitude')
    title(sprintf('%s',nameregion(i)))
end

SCATleadVSregionalmag2=figure(27)

stitle=suptitle('Scatter plot of lead time versus regional magnitude');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=5:8
    subplot(2,2,i-4)
    yesmag_region=yesmag(find(strcmp(yesregion,nameregion(i))));
    yeslead_region=yeslead(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesmag_region,yeslead_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesmag_region)>1
        [FIT,GOF]=fit(yesmag_region',yeslead_region','poly1');
        [Rval,Pval]=corrcoef(yesmag_region,yeslead_region);
        fitline=FIT.p1*yesmag_region+FIT.p2;
        plot(yesmag_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesmag),max(yeslead),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([min(yeslead) max(yeslead)])
    xlim([min(yesmag)-0.5 max(yesmag)+0.5])
    ylabel('Lead time (days)');
    xlabel('Magnitude')
    title(sprintf('%s',nameregion(i)))
end

% set(gcf, 'Position', get(0, 'Screensize'));
%% #28 & #29: Scatter plot of lead time versus regional depth

nameregion=unique(yesregion);
colpal=jet(numel(nameregion));

SCATleadVSregionaldep1=figure(28)

stitle=suptitle('Scatter plot of lead time versus regional depth');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=1:4
    subplot(2,2,i)
    yesdep_region=yesdep(find(strcmp(yesregion,nameregion(i))));
    yeslead_region=yeslead(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesdep_region,yeslead_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesdep_region)>1
        [FIT,GOF]=fit(yesdep_region',yeslead_region','poly1');
        [Rval,Pval]=corrcoef(yesdep_region,yeslead_region);
        fitline=FIT.p1*yesdep_region+FIT.p2;
        plot(yesdep_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesdep),max(yeslead),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([min(yeslead) max(yeslead)])
    xlim([min(yesdep)-10 max(yesdep)+10])
    ylabel('Lead time (days)');
    xlabel('Depth (km)')
    title(sprintf('%s',nameregion(i)))
end

SCATleadVSregionaldep2=figure(29)

stitle=suptitle('Scatter plot of lead time versus regional depth');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=5:8
    subplot(2,2,i-4)
    yesdep_region=yesdep(find(strcmp(yesregion,nameregion(i))));
    yeslead_region=yeslead(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesdep_region,yeslead_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesdep_region)>1
        [FIT,GOF]=fit(yesdep_region',yeslead_region','poly1');
        [Rval,Pval]=corrcoef(yesdep_region,yeslead_region);
        fitline=FIT.p1*yesdep_region+FIT.p2;
        plot(yesdep_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesdep),max(yeslead),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([min(yeslead) max(yeslead)])
    xlim([min(yesdep)-10 max(yesdep)+10])
    ylabel('Lead time (days)');
    xlabel('Depth (km)')
    title(sprintf('%s',nameregion(i)))
end

% set(gcf, 'Position', get(0, 'Screensize'));
%% #30 & #31: Scatter plot of lead time versus regional distance

nameregion=unique(yesregion);
colpal=jet(numel(nameregion));

SCATleadVSregionaldis1=figure(30)

stitle=suptitle('Scatter plot of lead time versus regional distance');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=1:4
    subplot(2,2,i)
    yesdis_region=yesdis(find(strcmp(yesregion,nameregion(i))));
    yeslead_region=yeslead(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesdis_region,yeslead_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesdis_region)>1
        [FIT,GOF]=fit(yesdis_region',yeslead_region','poly1');
        [Rval,Pval]=corrcoef(yesdis_region,yeslead_region);
        fitline=FIT.p1*yesdis_region+FIT.p2;
        plot(yesdis_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesdis),max(yeslead),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([min(yeslead) max(yeslead)])
    xlim([min(yesdis)-10 max(yesdis)+10])
    ylabel('Lead time (days)');
    xlabel('Distance (km)')
    title(sprintf('%s',nameregion(i)))
end

SCATleadVSregionaldis2=figure(31)

stitle=suptitle('Scatter plot of lead time versus regional distance');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=5:8
    subplot(2,2,i-4)
    yesdis_region=yesdis(find(strcmp(yesregion,nameregion(i))));
    yeslead_region=yeslead(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesdis_region,yeslead_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesdis_region)>1
        [FIT,GOF]=fit(yesdis_region',yeslead_region','poly1');
        [Rval,Pval]=corrcoef(yesdis_region,yeslead_region);
        fitline=FIT.p1*yesdis_region+FIT.p2;
        plot(yesdis_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesdis),max(yeslead),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([min(yeslead) max(yeslead)])
    xlim([min(yesdis)-10 max(yesdis)+10])
    ylabel('Lead time (days)');
    xlabel('Distance (km)')
    title(sprintf('%s',nameregion(i)))
end

% set(gcf, 'Position', get(0, 'Screensize'));
%% #32 Scatter plot of amplitude of adjusted Z/G (minus mu+2sig) versus magnitude, depth, distance and KLS
SCATadjzg1VSmagdepdisKLS=figure(32);

stitle=suptitle('Scatter plot of amplitude of adjusted Z/G (minus) versus magnitude, depth, distance and K_{LS}');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

subplot(2,2,1)
[FITmagZG,GOFmagZG]=fit(yesmag',yesadjZG_1','poly1');
magZG_fit=FITmagZG.p1*yesmag+FITmagZG.p2;
scatter(yesmag,yesadjZG_1,'bo','HandleVisibility','off');
hold on
plot(yesmag,magZG_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFmagZG.rsquare,GOFmagZG.rmse));
text(max(yesmag),max(yesadjZG_1),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Z/G');
xlabel('Magnitude')

subplot(2,2,2)
[FITdepZG,GOFdepZG]=fit(yesdep',yesadjZG_1','poly1');
depZG_fit=FITdepZG.p1*yesdep+FITdepZG.p2;
scatter(yesdep,yesadjZG_1,'bo','HandleVisibility','off');
hold on
plot(yesdep,depZG_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFdepZG.rsquare,GOFdepZG.rmse));
text(max(yesdep),max(yesadjZG_1),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Z/G');
xlabel('Depth (km)')

subplot(2,2,3)
[FITdisZG,GOFdisZG]=fit(yesdis',yesadjZG_1','poly1');
disZG_fit=FITdisZG.p1*yesdis+FITdisZG.p2;
scatter(yesdis,yesadjZG_1,'bo','HandleVisibility','off');
hold on
plot(yesdis,disZG_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFdisZG.rsquare,GOFdisZG.rmse));
text(max(yesdis),max(yesadjZG_1),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Z/G');
xlabel('Distance (km)')

subplot(2,2,4)
[FITKLSZG,GOFKLSZG]=fit(yesKLS',yesadjZG_1','poly1');
KLSZG_fit=FITKLSZG.p1*yesKLS+FITKLSZG.p2;
scatter(yesKLS,yesadjZG_1,'bo','HandleVisibility','off');
hold on
plot(yesKLS,KLSZG_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFKLSZG.rsquare,GOFKLSZG.rmse));
text(max(yesKLS),max(yesadjZG_1),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Z/G');
xlabel('K_{LS}')

% set(gcf, 'Position', get(0, 'Screensize'));
%% #33 Scatter plot of amplitude of adjusted Z/G (divide mu+2sig) versus magnitude, depth, distance and KLS
SCATadjzg2VSmagdepdisKLS=figure(33);

stitle=suptitle('Scatter plot of amplitude of adjusted Z/G (divide) versus magnitude, depth, distance and K_{LS}');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

subplot(2,2,1)
[FITmagZG,GOFmagZG]=fit(yesmag',yesadjZG_2','poly1');
magZG_fit=FITmagZG.p1*yesmag+FITmagZG.p2;
scatter(yesmag,yesadjZG_2,'bo','HandleVisibility','off');
hold on
plot(yesmag,magZG_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFmagZG.rsquare,GOFmagZG.rmse));
text(max(yesmag),max(yesadjZG_2),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Z/G');
xlabel('Magnitude')

subplot(2,2,2)
[FITdepZG,GOFdepZG]=fit(yesdep',yesadjZG_2','poly1');
depZG_fit=FITdepZG.p1*yesdep+FITdepZG.p2;
scatter(yesdep,yesadjZG_2,'bo','HandleVisibility','off');
hold on
plot(yesdep,depZG_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFdepZG.rsquare,GOFdepZG.rmse));
text(max(yesdep),max(yesadjZG_2),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Z/G');
xlabel('Depth (km)')

subplot(2,2,3)
[FITdisZG,GOFdisZG]=fit(yesdis',yesadjZG_2','poly1');
disZG_fit=FITdisZG.p1*yesdis+FITdisZG.p2;
scatter(yesdis,yesadjZG_2,'bo','HandleVisibility','off');
hold on
plot(yesdis,disZG_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFdisZG.rsquare,GOFdisZG.rmse));
text(max(yesdis),max(yesadjZG_2),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Z/G');
xlabel('Distance (km)')

subplot(2,2,4)
[FITKLSZG,GOFKLSZG]=fit(yesKLS',yesadjZG_2','poly1');
KLSZG_fit=FITKLSZG.p1*yesKLS+FITKLSZG.p2;
scatter(yesKLS,yesadjZG_2,'bo','HandleVisibility','off');
hold on
plot(yesKLS,KLSZG_fit)
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f',GOFKLSZG.rsquare,GOFKLSZG.rmse));
text(max(yesKLS),max(yesadjZG_2),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
hold off
legend('Linear fitting','location','northwest');
ylabel('Z/G');
xlabel('K_{LS}')

% set(gcf, 'Position', get(0, 'Screensize'));
%% #34 & #35: Scatter plot of amplitude of adjusted Z/G (divide mu+2sig) versus regional magnitude

nameregion=unique(yesregion);
colpal=jet(numel(nameregion));

SCATadjzg2VSregionalmag1=figure(34)

stitle=suptitle('Scatter plot of amplitude of adjusted Z/G (divide) versus regional magnitude');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=1:4
    subplot(2,2,i)
    yesmag_region=yesmag(find(strcmp(yesregion,nameregion(i))));
    yesadjZG_2_region=yesadjZG_2(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesmag_region,yesadjZG_2_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesmag_region)>1
        [FIT,GOF]=fit(yesmag_region',yesadjZG_2_region','poly1');
        [Rval,Pval]=corrcoef(yesmag_region,yesadjZG_2_region);
        fitline=FIT.p1*yesmag_region+FIT.p2;
        plot(yesmag_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesmag),max(yesadjZG_2),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([min(yesadjZG_2)-1 max(yesadjZG_2)+1])
    xlim([min(yesmag)-0.5 max(yesmag)+0.5])
    ylabel('Z/G');
    xlabel('Magnitude')
    title(sprintf('%s',nameregion(i)))
end

SCATadjzg2VSregionalmag2=figure(35)

stitle=suptitle('Scatter plot of amplitude of adjusted Z/G (divide) versus regional magnitude');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=5:8
    subplot(2,2,i-4)
    yesmag_region=yesmag(find(strcmp(yesregion,nameregion(i))));
    yesadjZG_2_region=yesadjZG_2(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesmag_region,yesadjZG_2_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesmag_region)>1
        [FIT,GOF]=fit(yesmag_region',yesadjZG_2_region','poly1');
        [Rval,Pval]=corrcoef(yesmag_region,yesadjZG_2_region);
        fitline=FIT.p1*yesmag_region+FIT.p2;
        plot(yesmag_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesmag),max(yesadjZG_2),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([min(yesadjZG_2)-1 max(yesadjZG_2)+1])
    xlim([min(yesmag)-0.5 max(yesmag)+0.5])
    ylabel('Z/G');
    xlabel('Magnitude')
    title(sprintf('%s',nameregion(i)))
end

% set(gcf, 'Position', get(0, 'Screensize'));
%% #36 & #37: Scatter plot of amplitude of adjusted Z/G (minus mu+2sig) versus regional magnitude

nameregion=unique(yesregion);
colpal=jet(numel(nameregion));

SCATadjzg1VSregionalmag1=figure(36)

stitle=suptitle('Scatter plot of amplitude of adjusted Z/G (minus) versus regional magnitude');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=1:4
    subplot(2,2,i)
    yesmag_region=yesmag(find(strcmp(yesregion,nameregion(i))));
    yesadjZG_1_region=yesadjZG_1(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesmag_region,yesadjZG_1_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesmag_region)>1
        [FIT,GOF]=fit(yesmag_region',yesadjZG_1_region','poly1');
        [Rval,Pval]=corrcoef(yesmag_region,yesadjZG_1_region);
        fitline=FIT.p1*yesmag_region+FIT.p2;
        plot(yesmag_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesmag),max(yesadjZG_1),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([min(yesadjZG_1)-1 max(yesadjZG_1)+1])
    xlim([min(yesmag)-0.5 max(yesmag)+0.5])
    ylabel('Z/G');
    xlabel('Magnitude')
    title(sprintf('%s',nameregion(i)))
end

SCATadjzg1VSregionalmag2=figure(37)

stitle=suptitle('Scatter plot of amplitude of adjusted Z/G (minus) versus regional magnitude');
pos_stitle=get(stitle,'Position');
set(stitle,'Position',[pos_stitle(1) pos_stitle(2)+0.02 pos_stitle(3)]);

for i=5:8
    subplot(2,2,i-4)
    yesmag_region=yesmag(find(strcmp(yesregion,nameregion(i))));
    yesadjZG_1_region=yesadjZG_1(find(strcmp(yesregion,nameregion(i))));
    hold on
    plot(yesmag_region,yesadjZG_1_region,'o','Color',colpal(i,:),'MarkerFaceColor',colpal(i,:),'HandleVisibility','off');
    if  numel(yesmag_region)>1
        [FIT,GOF]=fit(yesmag_region',yesadjZG_1_region','poly1');
        [Rval,Pval]=corrcoef(yesmag_region,yesadjZG_1_region);
        fitline=FIT.p1*yesmag_region+FIT.p2;
        plot(yesmag_region,fitline,'r')
        lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
        text(max(yesmag),max(yesadjZG_1),lbltxt,'HorizontalAlignment','right','VerticalAlignment','top');
        legend('Linear fitting','Location','northwest')
    end
    hold off
    ylim([min(yesadjZG_1)-1 max(yesadjZG_1)+1])
    xlim([min(yesmag)-0.5 max(yesmag)+0.5])
    ylabel('Z/G');
    xlabel('Magnitude')
    title(sprintf('%s',nameregion(i)))
end

% set(gcf, 'Position', get(0, 'Screensize'));
%% Saving figures

today=char(datetime('today','Format','dd-MM-yyyy'));
path='D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Statistik\Figures\';
completepath=strcat(path,today);
mkdir(completepath);

if exist('stitle')==1
    figname=string(stitle.String);
elseif exist('stitle')==0
    figname=get(gca,'Title');
    figname=string(figname.String);
end
figname=regexprep(figname,'/','');

numfig=numel(findobj('type','figure'));
for i=1:numfig
    set(gcf, 'Position', get(0, 'Screensize'));
    if numfig>1
        fignameext=strcat(figname,char(sprintf('-%d',i)));
        saveas(gcf,fullfile(completepath,fignameext),'tiff');
    else
        saveas(gcf,fullfile(completepath,figname),'tiff');
    end
    close gcf
end

clearvars stitle figname
