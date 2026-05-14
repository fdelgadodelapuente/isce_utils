clear all
close all
clc

filename='filt_220711-230612_8rlks_16alks_msk.unw.geo';
lambda='alos2';
zone='4V';
scaleval=0;
meanvel='no';
losfilename='220711-230612_8rlks_16alks.los.geo';
head_dir='new';
demf='crop.dem';
conncomp='filt_220905-230807_5rlks_28alks_msk.unw.conncomp.geo';
cor_file='220905-230807_5rlks_28alks.cor.geo';
cor_thresh=0;
coords.x1=NaN;

[ifg]=load_isce(filename,coords,lambda,zone,scaleval,...
    meanvel,losfilename,head_dir,demf);
    %%%%meanvel,losfilename,head_dir,demf,cor_file,cor_thresh,conncomp,4);
ulos=ifg.data;
X=ifg.X;
Y=ifg.Y;
x=ifg.lon;
y=ifg.lat;
z=ifg.z;
ny=ifg.ny;
nx=ifg.nx;

%
% ulos(1:200,900:nx)=NaN;

fwr=1; %scale factor for wrapping scale
phs_wr=ulos*4*pi/(ifg.lambda);
phs_wr=wrapToPi(phs_wr/fwr);
%phs_wr=wrapToPi(phs_wr);
phs_wr = (phs_wr/(2*pi) +1/2)*ifg.lambda*fwr/2;

figure;
    subplot(1,2,1)
    pcolor(ulos(1:10:end,1:10:end));shading flat;colormap jet;axis equal
    subplot(1,2,2)
	pcolor(phs_wr(1:10:end,1:10:end));shading flat;colormap jet;colorbar;axis equal

%%---- Ramp removal
%mask def sigmnal for ramp removal
ulos_ramp=ulos;
tmp=ulos;
x1b=900;
x2b=1300;
y1b=500;
y2b=800;
tmp(y1b:y2b,x1b:x2b)=NaN;

%define G matrix
xr=(X(:)-min(X(:)));
yr=(Y(:)-min(Y(:)));
zr=z(:)-min(z(:));
dcr=ones(size(Y(:)));
xryr=xr.*yr ;
%%%G=[ xr yr dcr zr]; %with topo corr, this particular data set can't be corrected with this
G=[ xr yr dcr];
G2=G;

%remove ramp from far-field data
phsr=tmp(:); 
ibs=isnan(G(:,1)) | isnan(phsr);
G(ibs,:)=[];
phsr(ibs)=[];
m=lsqlin(G,phsr);
ramp=G2*m;
ramp=reshape(ramp,ny,nx);
ramp(isnan(ulos))=NaN;
ulos=ulos-ramp; %remove ramp and save this
%%---- Ramp removal

indramp=665;
ulos_ax = [-0.35 0.75];

figure
    subplot(3,2,1)
        pcolor(x(1:5:end,1:5:end),y(1:5:end,1:5:end),ulos_ramp(1:5:end,1:5:end))
        hold on
        plot([x(1) x(end)],[y(indramp) y(indramp)],'k','LineWidth',2)
        shading flat
        axis equal
        colorbar
        %caxis([10 60])
        c=caxis;
        colormap jet
    subplot(3,2,2)
        hold on
        pcolor(x(1:5:end,1:5:end),y(1:5:end,1:5:end),tmp(1:5:end,1:5:end))
        plot([x(1) x(end)],[y(indramp) y(indramp)],'k','LineWidth',4)
        shading flat
        axis equal
        colorbar
        caxis(c)
        colormap jet
    subplot(3,2,3)
        hold on
        pcolor(x(1:5:end,1:5:end),y(1:5:end,1:5:end),ramp(1:5:end,1:5:end))
        plot([x(1) x(end)],[y(indramp) y(indramp)],'r','LineWidth',4)
        shading flat
        axis equal
        colorbar
        caxis(c)
        colormap jet
    subplot(3,2,4)
        hold on
        pcolor(x(1:5:end,1:5:end),y(1:5:end,1:5:end),ulos(1:5:end,1:5:end))
        plot([x(1) x(end)],[y(indramp) y(indramp)],'b','LineWidth',4)
        shading flat
        axis equal
        colorbar
        caxis(ulos_ax)
        colormap jet
    subplot(3,2,5:6)
        hold on
        plot(tmp(indramp,:),'k.')
        plot(ramp(indramp,:),'r-','LineWidth',4)
        plot(ulos(indramp,:),'b.')
        ylim(ulos_ax)
             
figure
        pcolor(x(1:10:end,1:10:end),y(1:10:end,1:10:end),ulos(1:10:end,1:10:end))
        shading flat
        axis equal
        colorbar
        caxis(ulos_ax)
        colormap jet

xmin=1; xmax=nx; ymin=1; ymax=ny;
%boxvar=phs(ymin:ymax,xmin:xmax);
%figure;pcolor(boxvar);shading flat;colormap jet;colorbar
%sqrt(nanvar(boxvar(:)))

ulos2wr = ulos;
phs_wr2 = wrapToPi(ulos2wr*4*pi/ifg.lambda);
phs_wr2 = (phs_wr2/(2*pi) +1/2)*ifg.lambda/2;



ramp_wr = wrapToPi(ramp*4*pi/ifg.lambda);
ramp_wr = (ramp_wr/(2*pi) +1/2)*ifg.lambda/2;

figure;pcolor(phs_wr2(1:3:end,1:3:end));shading flat;colormap jet; colorbar

ulos2export = ulos;
ulos2export (isnan(ulos2export )==1)=0;
sAs = single([flipud(ifg.mag) flipud(ulos2export*4*pi/ifg.lambda)]');
fid = fopen('filt_220910-240810_5rlks_28alks_msk_flat.unw.geo','wb');
fwrite(fid,sAs,'real*4');
fclose(fid);
             

ulos(isnan(ulos))=-100;
ulos(z<17)=NaN;
phs_wr(isnan(phs_wr))=-100;
phs_wr2(isnan(phs_wr2))=-100;
ramp_wr(isnan(ramp_wr))=-100;

grdwrite2(ifg.lon_vect(xmin:2:xmax),ifg.lat_vect(ymin:2:ymax),ulos(ymin:2:ymax,xmin:2:xmax)*100,strcat(filename,'.grd'));
%%grdwrite2(ifg.lon_vect(xmin:2:xmax),ifg.lat_vect(ymin:2:ymax),phs_wr(ymin:2:ymax,xmin:2:xmax)*100,strcat(filename,'wr.grd'));
%%grdwrite2(ifg.lon_vect(xmin:2:xmax),ifg.lat_vect(ymin:2:ymax),phs_wr2(ymin:2:ymax,xmin:2:xmax)*100,strcat(filename,'wr_flat.grd'));
%%grdwrite2(ifg.lon_vect(xmin:2:xmax),ifg.lat_vect(ymin:2:ymax),ramp_wr(ymin:2:ymax,xmin:2:xmax)*100,'ramp_wr.grd');
grdwrite2(ifg.lon_vect(xmin:2:xmax),ifg.lat_vect(ymin:2:ymax),z(ymin:2:ymax,xmin:2:xmax),'dem.grd');


 