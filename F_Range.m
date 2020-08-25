X=F_range(:,1);
Y=(F_range(:,1)+F_range(:,2))/2;
Z=F_range(:,2);
MAG=F_range(:,3);
DEPTH=F_range(:,4);
DIST=F_range(:,5);


j=1;
for i=1:26
   if F_range(i,1)==F_range(i,2)
       DIST2(j,1)=DIST(i,1);
       MAG2(j,1)=MAG(i,1);
       DEPTH2(j,1)=DEPTH(i,1);
       Y2(j,1)=F_range(i,1);
       j=j+1;
   end
end

j=1;
for i=1:26
   if DIST(i,1)>200
       DIST2(j,1)=DIST(i,1);
       Y2(j,1)=Y(i,1);
       j=j+1;
   end
end


figure(1)
for i=1:length(F_range)
    errorbar(F_range(i,3),(F_range(i,1)+F_range(i,2))/2,(F_range(i,2)-F_range(i,1))/2,'vertical','LineStyle','none','LineWidth',2,'Color',rand(1,3))
    hold on
end
xlim([3 10])
ylim([0 max(F_range(:,2))+0.01])
hold off
title('ULF frequency range against magnitude');
xlabel('Magnitude')
ylabel('ULF frequency range (Hz)')


figure(2)
for i=1:length(F_range)
    errorbar(F_range(i,4),(F_range(i,1)+F_range(i,2))/2,(F_range(i,2)-F_range(i,1))/2,'vertical','LineStyle','none','LineWidth',2,'Color',rand(1,3))
    hold on
end
ylim([0 max(F_range(:,2))+0.01])
hold off
title('ULF frequency range against depth');
xlabel('Depth (km)')
ylabel('ULF frequency range (Hz)')

figure(3)
for i=1:length(F_range)
    errorbar(F_range(i,5),(F_range(i,1)+F_range(i,2))/2,(F_range(i,2)-F_range(i,1))/2,'vertical','LineStyle','none','LineWidth',2,'Color',rand(1,3))
    hold on
end
ylim([0 max(F_range(:,2))+0.01])
hold off
title('ULF frequency range against epicentral distance');
xlabel('Distance (km)')
ylabel('ULF frequency range (Hz)')



j=0;
for i=1:length(F_range)
    if F_range(i,1)==F_range(i,2)
        j=j+1;
        F_del(j,1)=round(F_range(i,1),2);
        
    else
        x=round(F_range(i,1),2);
        z=round(F_range(i,2),2);
        for l=x:0.01:z
            j=j+1;
            F_del(j,1)=l;
        end
    end
end
F_del(F_del==0)=[];



for j=0.01:0.01:max(F_del)

        F_val(i,1)=j;
        F_val(i,2)=(sum(F_del==j)/numel(F_del))*100;
        j=j+0.01;
    
end
figure(4)
bar(F_val(:,1),F_val(:,2))
ylim([0 45]);
xlabel('Frequency (Hz)');
ylabel('Percentage of success (%)');
title('Percentage of success against frequency')
