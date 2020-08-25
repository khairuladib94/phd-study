%% Purpose
%To test statistical relationship of generated earthquake precursor data,
%interactively
%% Read excel
xlsfilename='D:\Study\Record of appearance of earthquakes geomagnetic precursor.xlsx';
xlssheetname='Adib_14g';
[readxcel1,readxcel2,readxcel]=xlsread(xlsfilename,xlssheetname);
today=char(datetime('today','Format','dd-MM-yyyy'));
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
%%
clc
close all

listvar=who('yes*');
listvar1=listvar(3);
disp('List of all available variables: ');
disp(listvar);
xinput=input('Input variable for x-axis: ','s');
yinput=input('Input variable for y-axis: ','s');
eval(['xdata=' xinput;]); eval(['ydata=' yinput;]);
xinput=upper(xinput(4:end)); yinput=upper(yinput(4:end));

clc
rem_unano=input('Remove unanomalous precursor? [y/n] ','s');
if rem_unano=='y'
    xdata(find(yesadjZG_2<1.0))=NaN;
    ydata(find(yesadjZG_2<1.0))=NaN;
end

clc
disp(sprintf('Categorization of data: \nNone (non)\nRegion (reg)\nCoordinates (coo)\nStation (stn)\nAnomalous (ano)\nYear (yer)'));
catdata=string(input('Input categorisation: ','s'));
listregion=unique(yesregion); liststn=unique(yesstn);

switch catdata
    case "non"
        listcat="All datapoints will be included.";
    case "reg"
        listcat=listregion;
        varcomp=yesregion;
    case "coo"
        listcat="Latitude: -90 to 90, longitude: -180 to 180";
    case "ano"
        listcat="Anomalous (1), not anomalous (0)";
    case "stn"
        listcat=liststn;
        varcomp=yesstn;
    case "yer"
        listcat=2005:2018;
        varcomp=string(yesyear);
end

clc
disp('Possible categorisation: ')
if catdata=="reg" || catdata=="stn"
    for i=1:numel(listcat)
        disp(sprintf('%d. %s',i,listcat(i)));
    end
else
    disp(listcat);
end

switch catdata
    case "reg"
        incldata=input('Input included categories. Input number(s) in [reg1,reg2,reg3...] format or 999 to include all: ');
    case "coo"
        incldata=input('Input included categories. Input number(s) in [lat1,lat2,lon1,lon2] format: ');
    case "ano"
        incldata=input('Input included categories. Input number(s) in [ano1,ano2] format: ');
    case "stn"
        incldata=input('Input included categories. Input number(s) in [stn1,stn2,stn3...] format or 999 to include all: ');
    case "yer"
        incldata=input('Input included categories. Input number(s) in [yer1,yer2,yer3...] format or 999 to include all: ');
end

if catdata=="reg"
    if incldata==999
        strindata=listregion;
    else
        strindata=listregion(incldata);
    end
elseif catdata=="coo"
    if incldata==999
        indata=find(yeslat);
    else
        indata=find(yeslat>=incldata(1) & yeslat<=incldata(2) & yeslon>=incldata(3) & yeslon<=incldata(4));
    end
elseif catdata=="stn"
    if incldata==999
        strindata=liststn;
    else
        strindata=liststn(incldata);
    end
elseif catdata=="yer"
    if incldata==999
        strindata=string(2005:2018);
    else
        strindata=string(incldata);
    end
end

fig=figure;
if catdata=="non"
    hold on
    plot(xdata,ydata,'bo','LineWidth',1.3);
elseif catdata=="reg" || catdata=="stn" || catdata=="yer"
    for i=1:numel(strindata)
        colpal=jet(numel(strindata));
        hold on
        plot(xdata(find(strcmp(varcomp,strindata(i)))),ydata(find(strcmp(varcomp,strindata(i)))),'o','Color',colpal(i,:),'LineWidth',1.3);
        if isempty(find(strcmp(varcomp,strindata(i)))) plot(NaN,NaN,'o','Color',colpal(i,:),'LineWidth',1.3); end
    end
    legend(strindata,'Location','northwestoutside');
elseif catdata=="coo"
    xxdata=xdata(indata); yydata=ydata(indata);
    plot(xxdata,yydata,'bo','LineWidth',1.3);
elseif catdata=="ano"
    colpal=jet(2);
    hold on
    if any(incldata==1) || incldata==999 
        plot(xdata(find(yesadjZG_2>1.0)),ydata(find(yesadjZG_2>1.0)),'o','Color',colpal(1,:),'LineWidth',1.3);
        if isempty(find(yesadjZG_2>1.0)) plot(NaN,NaN,'o','Color',colpal(1,:),'LineWidth',1.3); end
    end
    if any(incldata==0) || incldata==999 
        plot(xdata(find(yesadjZG_2<=1.0)),ydata(find(yesadjZG_2<=1.0)),'o','Color',colpal(2,:),'LineWidth',1.3);
        if isempty(find(yesadjZG_2<=1.0)) plot(NaN,NaN,'o','Color',colpal(2,:),'LineWidth',1.3); end
    end
    legend({'Anomalous','Not anomalous'},'Location','northwestoutside'); 
end
xlabel(xinput); ylabel(yinput);
figtitle=title(sprintf('Scatter plot of %s against %s',string(xinput),string(yinput)));
hold off

hold on
figobj=findobj(gca,'type','line');
xpts=[]; ypts=[];
for i=1:numel(figobj)
    xpt=figobj(i).XData; ypt=figobj(i).YData;
    xpts=horzcat(xpts,xpt); ypts=horzcat(ypts,ypt);
end
xpts(isnan(xpts))=[]; ypts(isnan(ypts))=[];
[FIT,GOF]=fit(xpts',ypts','poly1');
[Rval,Pval]=corrcoef(xpts,ypts);
fitline=FIT.p1*xpts+FIT.p2;
plot(xpts,fitline,'r','HandleVisibility','off');
lbltxt=string(sprintf('R^{2}=%.4f\nRMSE=%.4f\np-value=%.4f',GOF.rsquare,GOF.rmse,Pval(2,1)));
xlim=get(gca,'XLim'); ylim=get(gca,'YLim');
text('String',lbltxt,'HorizontalAlignment','right','VerticalAlignment','top','Units','normalized','EdgeColor','k','BackgroundColor','w','Position',[-0.048 0.085 0]);
lineeqn=string(sprintf('y=%.3fx%+.3f',FIT.p1,FIT.p2));
text('String',lineeqn,'HorizontalAlignment','right','VerticalAlignment','bottom','Position',[mean(xpts) mean(fitline)],'Color','r');
if rem_unano=='y'
    text('String','Unanomalous data removed.','HorizontalAlignment','right','VerticalAlignment','top','Units','normalized','EdgeColor','k','BackgroundColor','w','Position',[-0.048 0.5 0]);
end
hold off


clc
savefigure=input('Save figure? [y/n] ','s');
if savefigure=='y'
    set(gcf, 'Position', get(0, 'Screensize'));
    path='D:\OneDrive\Belajar\Dr. Adib Yusof\Data dan Analisis\Percubaan\Statistik\Figures\';
    completepath=strcat(path,today);
    mkdir(completepath);
    figname=figtitle.String;
    
    if exist(strcat(completepath,'\',figname,'.tif'))~=2
        saveas(gcf,fullfile(completepath,figname),'tiff');
    else
        for i=1:100
            fignameext=strcat(figname,char(sprintf('-%d',i)));
            if exist(strcat(completepath,'\',fignameext,'.tif'))~=2
                saveas(gcf,fullfile(completepath,fignameext),'tif');
                break;
            end
        end
    end
    disp('Figure saved');
end
clc
disp('Analysis ends');
pause(0.5)
clc










