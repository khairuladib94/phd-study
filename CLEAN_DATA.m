
H_dn=H;
D_dn=D;
Z_dn=Z;
F_dn=F;
%%

% H_dn=medfilt1(H_dn,30,'omitnan','truncate');
% D_dn=medfilt1(D_dn,30,'omitnan','truncate');
Z_dn=medfilt1(Z_dn,30,'omitnan','truncate');
F_dn=medfilt1(F_dn,30,'omitnan','truncate');

%%

    
H_sig=std(H_dn,'omitnan');
H_mu=mean(H_dn,'omitnan');
D_sig=std(D_dn,'omitnan');
D_mu=mean(D_dn,'omitnan');
Z_sig=std(Z_dn,'omitnan');
Z_mu=mean(Z_dn,'omitnan');
F_sig=std(F_dn,'omitnan');
F_mu=mean(F_dn,'omitnan');

for j=1:length(H_dn)
    if H_dn(j)>H_mu+5*H_sig||H_dn(j)<H_mu-5*H_sig
%         H_dn(j)=NaN;
    end
    if D_dn(j)>D_mu+5*D_sig||D_dn(j)<D_mu-5*D_sig
%         D_dn(j)=NaN;
    end
    if Z_dn(j)>Z_mu+5*Z_sig||Z_dn(j)<Z_mu-5*Z_sig
        Z_dn(j)=NaN;
    end
    if F_dn(j)>F_mu+5*F_sig||F_dn(j)<F_mu-5*F_sig
        F_dn(j)=NaN;
    end
end

%%
figure
sp1=subplot(4,1,1);
plot(UT1m,H_dn)
title('H')

sp2=subplot(4,1,2);
plot(UT1m,D_dn)
title('D')

sp3=subplot(4,1,3);
plot(UT1m,Z_dn)
title('Z')

sp4=subplot(4,1,4);
plot(UT1m,F_dn)
title('F')

linkaxes([sp1,sp2,sp3,sp4],'x')











