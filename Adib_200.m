%% Purpose
% Observe Pz/Pg in 9 frequency bands, e.g. 0.01 - 0.02, 0.02 - 0.03 .. 0.09
% - 0.10 hz
% Estimate angle using polarization ellipse method as used by Schekotov
% The frequency band to bandpass filter horizontal components to calculate
% azimuthal angle can be chosen, e.g. 0.01 - 0.02, 0.02 - 0.03 .. 0.09 - 0.10 hz
% 2007, 2008, 2015 and Ohta 2013
% Altered from Adib_104
%% Load initial variable
clc; 
close all;
load VARIABLES_WORLD
%%
UTC = datetime(EQ_WORLD(:, 1:6));
[Lat, Lon, Depth, Mag]= deal(EQ_WORLD(:, 7), EQ_WORLD(:, 8), EQ_WORLD(:, 9), EQ_WORLD(:, 10));
EQTab = table(UTC, Lat, Lon, Depth, Mag);

LT_start=[22,23,00,01,02];
LT_end=  [02,03,04,05,06];
StudiedStn = string(StudiedStn);
ZGTab = table;
PeriodCount = 1;

for iStn = 1:numel(StudiedStn)
    stn_latlon=[stn_MAGDAS(StudiedStn(iStn) == string(stn_vec),2:3)];
    EQ_dis=deg2km(distance('gc',stn_latlon,[EQTab.Lat, EQTab.Lon]));
    EQDateNearStn = dateshift(EQTab.UTC(EQ_dis < 180), 'start', 'day');
    
    for iYear = 2007:2016
        disp(string(StudiedStn(iStn)) + string(iYear) + "...");
        filepathname=fullfile('E:\Study\MAGDAS data\',StudiedStn(iStn),'\',StudiedStn(iStn)) + string(iYear) + 'S.mat';
        clearvars H D Z F UT1m
        if exist(filepathname,'file'), load(filepathname); else, continue; end
        UT1m = (datetime(UT1m, 'ConvertFrom', 'datenum'))';
        [G_dn, Z_dn] = preProcess(H, D, Z);
        
        clearvars NonEQPeriods
        date_start=datetime(iYear, 1, 1);
        IsNextYear = false;
        i = 1;
        while ~ IsNextYear
            date_end=date_start + 60;
            IsEQInPeriod = date_start <= EQDateNearStn & date_end >= EQDateNearStn;
            if all(~IsEQInPeriod)
                NonEQPeriods(i, :) = [date_start, date_end];
                date_start = date_end + 1;
                i = i + 1;
            else
                date_start = EQDateNearStn(find(IsEQInPeriod, 1, 'last')) + 1;
            end
            IsNextYear = year(date_start + 60) > iYear;
        end
        if ~ exist('NonEQPeriods', 'var') || isempty(NonEQPeriods), continue; end
        
        for iPeriod = 1:height(NonEQPeriods)
            GappedObs = false;
            DataStart = find(UT1m == NonEQPeriods(iPeriod, 1));
            DataEnd = find(UT1m == NonEQPeriods(iPeriod, 2)+1-seconds(1) );
            Tmp = DataStart:DataEnd;
            [G_period, Z_period, UT1m_period] = deal(G_dn(Tmp), Z_dn(Tmp), UT1m(Tmp));
            [ap_daily, Dst_daily] = getIndex(index_geomag, UT1m_period);
            [ap_mean, Dst_mean] = deal(mean(ap_daily, 'omitnan'), mean(Dst_daily, 'omitnan'));
            for iLT = 1:5
                [LTStart, LTEnd] = deal(LT_start(iLT), LT_end(iLT));
                [Z_night, G_night]  = extractNighttime(G_period, Z_period, UT1m_period, LTStart, LTEnd, stn_latlon);
                
                start_date = NonEQPeriods(iPeriod, 1);
                end_date = NonEQPeriods(iPeriod, 2);
                ZG = PRA(Z_night, G_night);
                if sum(isnan(ZG(1,:))) > 15, GappedObs = true; end
                if GappedObs, break; end
                clearvars ZGAnom
                for iFreq = 1:9
                    ZGAnom(iFreq, :) =  ZG(iFreq, :) > mean(ZG(iFreq, :), 'omitnan') + 2*std(ZG(iFreq, :), 'omitnan') & ...
                        ZG(iFreq, :) > 1.5 & ap_mean <= 27 & Dst_mean >= -30;
                end
                ZGTab.Stn(PeriodCount:PeriodCount+8) = StudiedStn(iStn);
                ZGTab.StartDate(PeriodCount:PeriodCount+8) = NonEQPeriods(iPeriod, 1);
                ZGTab.EndDate(PeriodCount:PeriodCount+8) = NonEQPeriods(iPeriod, 2);
                ZGTab.LT(PeriodCount:PeriodCount+8) = iLT;
                ZGTab.Freq(PeriodCount:PeriodCount+8) = 1:9;
                for iFreq = 1:9
                    ZGTab.Anom(PeriodCount+iFreq-1) = {ZGAnom(iFreq, :)};
                    ZGTab.AnomCount(PeriodCount+iFreq-1) = sum(ZGAnom(iFreq, :));
                end
                PeriodCount = PeriodCount + 9;
            end
            
        end  
    end
end
%% Plotting
CumAnom = zeros(1,61);
for i = 1:height(ZGTab)
   if ZGTab.AnomCount(i) == 0, continue; end
   CumAnom = CumAnom + ZGTab.Anom{i, :};
end
PercentAnom = 100*CumAnom / height(ZGTab);
bar(PercentAnom, 'FaceColor', [0.80,0.80,0.80]);
title(['$n = ', num2str(height(ZGTab)), '$']);
ylabel('Percentage of anomalies (\%)'); xlabel('Day');
xlim([0.5, 61.5]);
xticks(1:5:61);
Ax = gca;
nicefigure;
%% Local functions
function [H_dn, G_dn] = preProcess(H, D, Z)
H_dn=medfilt1(H,1,'omitnan','truncate');
D_dn=medfilt1(D,1,'omitnan','truncate');
Z_dn=medfilt1(Z,1,'omitnan','truncate');

IsOutlier = @(X) abs(X) > mean(X, 'omitnan') + 7*std(X, 'omitnan');
H_dn(IsOutlier(H_dn))=NaN;
D_dn(IsOutlier(D_dn))=NaN;
Z_dn(IsOutlier(Z_dn))=NaN;
G_dn = sqrt(H_dn.^2 + D_dn.^2);
end

function [Z_night, G_night] = extractNighttime(G_period, Z_period, UT1m_period, LTStart, LTEnd, stn_latlon)

LT_length=LTEnd-LTStart;
if LT_length<0 LT_length=LT_length+24; end
t_zone=timezone(stn_latlon(2),'degrees');
t_start=LTStart+t_zone;
if t_start>24 t_start=t_start-24; end
if t_start<0 t_start=t_start+24;end
t_starts=t_start*3600;
t_int=LT_length*3600;

days_num = int8(days(range(UT1m_period)));
for i=1:days_num
    Z_night(:,i)=Z_period(t_starts+1:t_starts+t_int);
    G_night(:,i)=G_period(t_starts+1:t_starts+t_int);
    UT1m_night(:,i)=UT1m_period(t_starts+1:t_starts+t_int);
    t_starts=t_starts+86400;
end

%Gap-filling
for i=1:days_num
    Z_measure=regionprops(logical(isnan(Z_night(:,i))), 'Area');
    Z_length=[Z_measure.Area]; if isempty(Z_length) Z_length=1; end
    if max(Z_length)<0.01*size(Z_night,1) && any(isnan(Z_night(:,i)))
        Z_night(:,i)=fillgaps(Z_night(:,i));
    end
    
    G_measure=regionprops(logical(isnan(G_night(:,i))), 'Area');
    G_length=[G_measure.Area]; if isempty(G_length) G_length=1; end
    if max(G_length)<0.01*size(G_night,1) && any(isnan(G_night(:,i)))
        G_night(:,i)=fillgaps(G_night(:,i));
    end
end
end

function [ap_daily, Dst_daily] = getIndex(index_geomag, UT1m_period)
index_geomag_Idx = find(ismember( datetime(index_geomag(:,1), 1, index_geomag(:,2), index_geomag(:,3), 0, 0), ...
    UT1m_period(1):hours(1):dateshift(UT1m_period(end), 'start', 'hour')));
[ap, Dst]=deal(index_geomag(index_geomag_Idx,6), index_geomag(index_geomag_Idx,5));
days_num = int8(days(range(UT1m_period)));
ap_daily=reshape(ap,[],days_num);
Dst_daily=reshape(Dst,[],days_num);
end

function ZG = PRA(Z_night, G_night)
f_rangeint=[0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10];
days_num = 61;
alpha1=NaN(size(G_night,1),days_num);
alpha2=NaN(size(G_night,1),days_num);
for i=1:days_num
    
    % Skip if data contains nans
    if any(isnan(Z_night(:,i))) || any(isnan(G_night(:,i)))
        ZZ(:,i)=NaN(1,9);
        GG(:,i)=NaN(1,9);
        continue
    end
    
    % Power spectral density
    n_length=size(G_night,1);
    n_win=n_length;
    n_ovrlp=0;
    n_fft=n_win;
    fs=1;
    
    [Pzz,fff]=periodogram(Z_night(:,i),hamming(n_win),n_fft,fs);
    [Pgg,fff]=periodogram(G_night(:,i),hamming(n_win),n_fft,fs);
    
    bin_floor=find(round(fff,4)==f_rangeint(1)); bin_ceil=find(round(fff,4)==f_rangeint(end));
    
    % ULF extraction
    Pzz1=Pzz(bin_floor:bin_ceil,:);
    Pgg1=Pgg(bin_floor:bin_ceil,:);
    fff1=fff(bin_floor:bin_ceil,:);
    
    % ULF averaging
    for b=1:numel(f_rangeint) f_rangeidx(b)=find(round(fff1,4)==f_rangeint(b)); end
    
    for a=1:numel(f_rangeint)-1
        Pzz_band(:,a)=Pzz1(f_rangeidx(a):f_rangeidx(a+1)-1);
        Pgg_band(:,a)=Pgg1(f_rangeidx(a):f_rangeidx(a+1)-1);
        
        fff_band(:,a)=fff1(f_rangeidx(a):f_rangeidx(a+1)-1);
    end
    
    % Polarization ratio calculation
    ZZ(:,i)=real(mean(Pzz_band));
    GG(:,i)=real(mean(Pgg_band));
    
end
ZZ_norm=normalize(ZZ,2,'range',[1 3]);
GG_norm=normalize(GG,2,'range',[1 2]);
ZG=ZZ_norm./GG_norm;

end
