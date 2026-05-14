%decompose descending and ascending U_los into U_x and U_z 
%interferograms must be in the JPL RMG file format (stripmapApp.py,
%topsApp.py, alos2App.py), and geocoded to the same bounding box
%
%Francisco Delgado, Cornell University
%
%December 2015 for Cordon Caulle RS2/CSK
%Updated September 2021 for Cumbre Vieja S1
%Updated December 2022 for Mauna Loa S1/CSK

close all
clear all
clc

zone ='18G';
ind2plot=645;
axcords=[800 1400 400 800 ];

fn_asc='alos2/pairs/alos2_asc_2022_2024.geo';
wvl_asc='alos2';
losf_asc='alos2/pairs/191106_8rlks_16alks.los.geo';
fn_dsc='saocom_dsc/Igrams/saocom_dsc_2023_2024.geo';
wvl_dsc='saocom';
losf_dsc='saocom_dsc/Igrams/los.rdr.geo';

coords.x1 = NaN  ;
asc=load_isce(fn_asc,coords,wvl_asc,zone,0,'no',losf_asc,'new');
asc.data=[-asc.data];
asc_s1=asc.S(:,:,1);
asc_s2=asc.S(:,:,2);
asc_s3=asc.S(:,:,3);
asc_s1=asc_s1(:);
asc_s2=asc_s2(:);
asc_s3=asc_s3(:);
ascv=asc.data(:);
figure;pcolor(asc.data(1:10:end,1:10:end));shading flat;caxis([-0.05 0.2]);colormap jet
boxx = asc.data(440:450,1270:1280);
%ascv = ascv - nanmean(boxx(:));

coords.x1 = NaN;
dsc=load_isce(fn_dsc,coords,wvl_dsc,zone,0,'no',losf_dsc,'new');
dsc.data=[-dsc.data];
dsc_s1=dsc.S(:,:,1);
dsc_s2=dsc.S(:,:,2);
dsc_s3=dsc.S(:,:,3);
dscv=dsc.data;

dscv=griddata(dsc.X,dsc.Y,dscv,asc.X,asc.Y);
dsc_s1=griddata(dsc.X,dsc.Y,dsc_s1,asc.X,asc.Y);
dsc_s3=griddata(dsc.X,dsc.Y,dsc_s3,asc.X,asc.Y);
dscv_resized=dscv;

figure;pcolor(dsc.data(1:1:end,1:1:end));shading flat;caxis([-0.05 0.2]);colormap jet
boxx = dscv(440:450,1270:1280);
%dscv = dscv - nanmean(boxx(:));

%%%figure; pcolor(asc.data); shading flat; axis equal
%%%figure; pcolor(dscv); shading flat; axis equal

dsc_s1=dsc_s1(:);
dsc_s2=dsc_s2(:);
dsc_s3=dsc_s3(:);
dscv=dscv(:);

%account for the variance in the different data sets
sig_a = 0.01; %asc variance
sig_d = 0.01; %dsc variance
covd=diag([sig_a sig_d]);
ch    = chol(covd);
Cdinv = inv(ch');

for i =1:length(ascv)
    if (isnan(ascv(i))==0 & isnan(dscv(i))==0) %select pixels with 2 measurements
        G=-[asc_s1(i) asc_s3(i); dsc_s1(i) dsc_s3(i)]; % (-) to account for the look instead of the heading vector
        d=[ascv(i); dscv(i)];
        Gw=Cdinv*G;
        dw=Cdinv*d;
        m=Gw\dw;
        %m=inv(G'*inv(covd)*G)*G'*inv(covd)*d;
        %C=diag(inv(G'*inv(covd)*G));
        ux(i)=m(1);
        uz(i)=m(2);
        rat(i)=ux(i)/uz(i);
    else
        ux(i)=NaN;
        uz(i)=NaN;
        rat(i)=NaN;
    end
end
        
ux=reshape(ux,[size(asc.data)]);
uz=reshape(uz,[size(asc.data)]);
rat=reshape(rat,[size(asc.data)]);

ux2=ux;
uz2=uz;
ux2(ux2==0)=1;
uz2(uz2==0)=1;

%get minimum
k=find(min(ux(:))==ux(:));
[min_j,min_i]=ind2sub(size(dsc.X),k);

ratio_max=max(abs(ux2(:)))/abs(max(uz2(:)));
disp(['ratio max = ',num2str(ratio_max)])
ratio_min=min((ux2(:)))/(max(uz2(:)));
disp(['ratio min = ',num2str(ratio_min)])
disp(['Uz max = ',num2str(max(uz(:)))]);

%find zone of max uplift
Uz=uz*0;
%%%%%%Uz(550:700,500:800)=uz(550:700,500:800);%window for the maximum
Uz = uz;
Uz_filt = medfilt2(Uz,[3 3]);
k=find(max(Uz_filt(:))==Uz_filt(:));
ascx=asc.X(:);
ascy=asc.Y(:);
disp(['Max uplift ',num2str(round([ascx(k) ascy(k)]))])
[l,l2]=my_utm2ll(ascx(k),ascy(k),1,'18G');
[j,i]=ind2sub(size(asc.X),k);
uz_max=[l l2];
disp([num2str([l l2])]);
save('uz_max_2024.txt','uz_max','-ascii')

xcord=asc.X(ind2plot,:)-axcords(1)*1e3;

px=[asc.X(ind2plot,1) asc.X(ind2plot,end)]/1e3;
py=[asc.Y(ind2plot,1) asc.Y(ind2plot,length(asc.Y(ind2plot,:)))]/1e3;

xprof=[200 400]*2;
yprof=[250 400]*2;
ux_perf=improfile(ux,xprof,yprof);
uz_perf=improfile(uz,xprof,yprof);

figure;
    subplot(3,2,1)
        colormap jet
        hold on;
        pcolor(asc.data*100);
        shading flat;
        plot(min_i,min_j,'k.','MarkerSize',20);
        plot(xprof,yprof)
        colorbar
        %axis(axcords)
        caxis([-20 20])
        axis equal
    subplot(3,2,2)
        colormap jet
        hold on;
        pcolor(dscv_resized*100);
        shading flat;
        plot(min_i,min_j,'k.','MarkerSize',20);
        plot(xprof,yprof)
        colorbar
        %axis(axcords)
        caxis([-20 20])
        axis equal
    subplot(3,2,3)
        colormap jet
        hold on;
        pcolor(ux*100);
        shading flat;
        plot(min_i,min_j,'k.','MarkerSize',20);
        plot(xprof,yprof)
        colorbar
        %axis(axcords)
        caxis([-20 20])
        axis equal
    subplot(3,2,4)
        colormap jet
        hold on;
        pcolor(uz*100);
        shading flat;
        plot(i,j,'k.','MarkerSize',20);
        plot(xprof,yprof)
        colorbar
        %axis(axcords)
        caxis([-20 20])
        axis equal
    subplot(3,2,5)
        hold on
        plot(ux_perf*100,'.k')
        %plot(uz_perf/max(uz_perf(:)),'.b')
    subplot(3,2,6)
        hold on
        %plot(ux_perf/max(uz_perf(:)),'.k')
         plot(uz_perf*100,'.b')

        
wvl=asc.lambda;
uz(isnan(uz)==1)=0;
sAs = single([flipud(0*uz) flipud(uz*4*pi/wvl)]');
fid = fopen('uz.geo','wb');
fwrite(fid,sAs,'real*4');
fclose(fid);   


ux(isnan(ux)==1)=0;
sAs = single([flipud(0*ux) flipud(ux*4*pi/wvl)]');
fid = fopen('ux.geo','wb');
fwrite(fid,sAs,'real*4');
fclose(fid);   

ux2ex=ux;
uz2ex=uz;
ux2ex(ux2ex==0)=-99;
uz2ex(uz2ex==0)=-99;
grdwrite2(asc.lon_vect,asc.lat_vect,ux2ex*100,'ux.grd')
grdwrite2(asc.lon_vect,asc.lat_vect,uz2ex*100,'uz.grd')