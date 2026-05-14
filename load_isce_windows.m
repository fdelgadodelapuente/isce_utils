function [datastruct]=load_isce(filename,coords,lambda,zone,scaleval,...
    meanvel,losfilename,head_dir,demf,cor_file,cor_thresh,conncomp,conncomp_ind)

%ISCE interferogram loader, 
%
%Francisco Delgado, Cornell University, Institut de Physique du Globe de
%Paris
%
%July 2017 original script, modified from Rowena Lohman's load_any_data for
%ROI_PAC
%June 2019, parse the coordinates directly from .vrt file
%January 2020, parse coherence file reading number of bands in the .vrt file
%
% [datastruct]=load_isce(filename,coords,lambda,zone,scaleval,...
%    meanvel,losfilename,head_dir,demf,cor_file,cor_thresh)
%
%filename       interferogram
%coords         data structure with coordinates, deprecated. use [] to read a .vrt file
%lambda         satellite: alos2, alos, csk, tsx, envisat, sentinel1, rs2
%zone           utm zone
%scaleval       zero padding scale factor for downsampling, 0 for anything else
%meanvel        'yes' for mean velocity from time series, 'no' for interferograms
%losfilename    line-of-sight file
%head_dir       heading convention, either 'old' or 'new' 
%demf           DEM file, to be used for ramp removal or for plotting
%cor_file       correlation file
%cor_thresh     mask pixels below this threshold
%
% 
%example
%[datastruct]=load_isce('filt_topophase.unw.geo',[],'tsx','19G',0,'no',...
%'los.rdr.geo','new','dem.crop','topophase.cor.geo',.1);

if(nargin<5)
    scaleval=0;
end
% if(nargin<8)
%     head_dir='new';
% else
%     head_dir='old';
% end
if isfinite(coords.x1)==0 %%parse coordinates from .vrt file
% %     cmd = strcat('grep rasterXSize ',{' '},filename,'.vrt');
% %     [b,tmp]=unix(cmd{1});
% %     %nx=str2double(tmp(26:29));
% %     %ny=str2double(tmp(45:47));
% %     tmp=str2double(regexp(tmp,'[\d.]+','match'));
% %     nx=tmp(1);
% %     ny=tmp(2);
% % 
% %     cmd = strcat('grep GeoTransform ',{' '},filename,'.vrt');
% %     [b,tmp]=unix(cmd{1});
% %     tmp=convertCharsToStrings(tmp);
% %     tmp = erase(tmp,'<GeoTransform>');
% %     tmp = erase(tmp,'</GeoTransform>');
% %     tmp = erase(tmp,',');
% %     tmp = split(tmp,[" "]);
% %     tmp = tmp(end-5:end);
% %     x1=str2num(tmp(1));
% %     dx=str2num(tmp(2));
% %     y2=str2num(tmp(4));
% %     dy=str2num(tmp(6));
else
    x1=coords.x1;
    y2=coords.y2;
    dx=coords.dx;
    nx=coords.nx;
    ny=coords.ny;
    dy=-dx;
end

%satellites wvls
wvl.csk       = 0.031228381041666666;
wvl.tsx       = 0.0310665781638;
wvl.envisat   = 0.0562356424;
wvl.rs2       = 0.055465772432999993;
wvl.sentinel  = 0.05546576;
wvl.alos      = 0.236057;
wvl.alos2     = 0.2424525;
wvl.dem       = 1; %for DEMs, no scaling

disp(['Opening ',filename,' satellite wvl ',lambda])
disp(['x1 ','nx ', 'dx ', 'y2 ', 'ny ', 'dy '])
disp([num2str(x1,5) num2str(nx) num2str(dx) num2str(y2) num2str(ny) num2str(dy), {' '}])

if strcmp(lambda,'alos2')==1;     lambda =wvl.alos2; end
if strcmp(lambda,'alos')==1;      lambda =wvl.alos;  end
if strcmp(lambda,'csk')==1;       lambda =wvl.csk;   end
if strcmp(lambda,'tsx')==1;       lambda =wvl.tsx;   end
if strcmp(lambda,'envisat')==1;   lambda =wvl.envisat; end
if strcmp(lambda,'sentinel')==1;  lambda =wvl.sentinel; end
if strcmp(lambda,'rs2')==1;       lambda =wvl.rs2; end
if strcmp(lambda,'dem')==1;       lambda =wvl.dem; end %to compensate for line 75
filename
h =fopen(filename,'r');
[F,count] = fread(h,2*nx*ny,'float32');
rmg = reshape(F,2*nx,ny);
phs = rmg(nx+1:nx*2,:);
mag = flipud(rmg(1:nx,:)');
phs = flipud(phs'); % need to switch rows and columns
phs(phs==0)=NaN;
%phs(mag<1)=NaN;
%phs(mag<1e3)=NaN; %mask edges with filtering stripes

if strcmp(meanvel,'yes')==1;
    data=-phs/100; %cm to m
else if strcmp(meanvel,'dem')==1;
	data=-phs;
else
    data  = -phs*lambda/(4*pi) ; %phase to meters
    end
end

if(scaleval)
    olddata = data;
    scale   = 2^ceil(log2(min([ny nx])/scaleval));
    extrax  = scale-mod(nx,scale);
    extray  = scale-mod(ny,scale);
    data    = blkdiag(olddata,zeros(extray,extrax));
    [ny,nx] = size(data);
    y2      = y2-extray*dy;
    disp(['padding X,Y,data,S with ' num2str(extrax) '/' num2str(extray) ' zeros']);
else
    extrax = 0;
    extray = 0;
    scale  = 0;
end

%get coordinates from upper left vertex
y1=y2+dy*(ny-1);
x2=x1+dx*(nx-1);

x=x1:dx:x2;
y=y1:-dy:y2;

[X,Y]=meshgrid(x,y);
lon=X;
lat=Y;

%convert to cartesian coordinates
%[X,Y]=my_utm2ll(X,Y,2,zone); %don't use default function, only the mappingtoolbox
    r_a    = 6378137.0;            %equatorial radius, WGS84
    r_pole = 6356752.3;            %polar radius, WGS84
    r_e2   = (1-r_pole^2/r_a^2);   %ellipsoid eccentricity squared
    r_e    = sqrt(r_e2);
    mstruct       = defaultm('utm');
    mstruct.geoid = [r_a r_e];
    mstruct.zone = zone;
    mstruct      = defaultm(mstruct);
    [X,Y]        = mfwdtran(mstruct,Y,X);

pixelsize=mean([sqrt((X(1)-X(2))^2+(Y(1)-Y(2))^2) sqrt((X(nx*ny-1)-X(nx*ny))^2+(Y(nx*ny-1)-Y(nx*ny))^2)]);
X=X+pixelsize/2;%We use a coord system with the upper left corner of each pixel but the dem coord in the dem.rsc file is usually the center.
Y=Y-pixelsize/2;    

%look for bad points
baddata         = mode(data(:));
badid           = find(data==baddata);
data(badid)     = NaN;
disp(['setting ' num2str(length(badid)) ' pts with phs=' num2str(baddata) ' to NaN']);

datastruct=struct('data',data,'mag',mag,'phs',phs,'X',X,'Y',Y,'pixelsize',pixelsize, ...
    'zone',zone,'lambda',lambda,'nx',nx,'ny',ny,'filename',filename, ...
    'scale',scale','extrax',extrax,'extray',extray,'lon',lon,'lat',lat,'lon_vect',x,'lat_vect',y);

if(nargin>=7)

        %Load LOS
        extrax = datastruct.extrax;
        extray = datastruct.extray;
        ox     = nx-extrax;
        oy     = ny-extray;

        fid          = fopen(losfilename,'r','native');
        [temp,count] = fread(fid,[ox*2,oy],'real*4');
        status       = fclose(fid);

        look = temp(1:ox,:); %first part of file
        head = temp(ox+1:ox*2,:); %second part of file
        head = flipud(head'); %this is right
        look = flipud(look');%this is also right

        squint  = 0; %Squint is already taken into account in the los.rdr file

        if strcmp(head_dir,'old')==1;
            head(head~=0) = -head(head~=0)+180; %clockwise for isce201407, CSK 2012-2015 Caulle - 2015 Villarrica
        else if strcmp(head_dir,'new')==1;
            head(head~=0) = -(head(head~=0)+180);  %counterclockwise for rs2App, alos2App and isce201604, 201609, 201704, zd_isce201506
        end
        end
        
        head = (head-squint).*pi/180;
        look    = look.*pi/180;

        id          = find(head==0);
        jd          = find(head~=0);
        head(id) = mean(head(jd));
        look(id)    = mean(look(jd));

        S1 = [sin(head).*sin(look)];
        S2 = [cos(head).*sin(look)];
        S3 = [ -cos(look)];

        S1      = blkdiag(S1,zeros(extray,extrax)); %zero padding for downsampling
        S2      = blkdiag(S2,zeros(extray,extrax)); %zero padding for downsampling
        S3      = blkdiag(S3,zeros(extray,extrax)); %zero padding for downsampling
        badid   = find(S1(:)==0);
        S1(badid) = S1(1); % set to average in load_los
        S2(badid) = S2(1);
        S3(badid) = S2(1);

        if strcmp(meanvel,'dem')==1;  %only vertical displacement
            S1=S1*0;
            S2=S2*0;
            S3=S3*0-1;
        end
        
        S(:,:,1)  = S1;
        S(:,:,2)  = S2;
        S(:,:,3)  = S3;
        
        disp('average angles')
        head2=head;         head2(head==0)=NaN;
        look2=look;         look2(head==0)=NaN;
        disp([nanmean(head2(:)*180/pi) nanmean(look2(:)*180/pi)])
        disp('average LOS')
        S1_2=S1;         S1_2(gradient(S1_2)==0)=NaN;
        S2_2=S2;         S2_2(gradient(S2_2)==0)=NaN;
        S3_2=S3;         S3_2(gradient(S3_2)==0)=NaN;
        disp([nanmean(S1_2(:)) nanmean(S2_2(:)) nanmean(S3_2(:))])

        datastruct.S=S;

        datastruct.head      = blkdiag(head,zeros(extray,extrax));
        datastruct.look      = blkdiag(look,zeros(extray,extrax));
end

if(nargin>=10)
    cor_file0=cor_file;
%     if strcmp(cor_file0(1),'2')==1
%         cor_file0=cor_file0(19:end); %remove prefix for TS
%     end
%     %cor_file=cor_file(19:end); %for caulle 2016-2018 paper??
%     %%%%if strcmp(cor_file(19:21),'phs')==1 %phsig.cor, only for zerodop old file format
%     cmd = strcat('grep band ',{' '},cor_file,'.vrt | wc -l');
%     [~,bands]=unix(cmd{1}); bands=str2num(bands);
%     %%%%if strcmp(cor_file0(1:3),'phs')==1 %phsig.cor, only for zerodop old file format
%     if bands==1 %phsig.cor, use modern .vrt files
        h =fopen(cor_file,'r');
        [F,count] = fread(h,nx*ny,'float32');
        %[F,count] = fread(h,nx*ny,'real*4');
        cor = reshape(F,nx,ny);
        cor = flipud(cor'); % need to switch rows and columns
%     else %topophase.cor
%         fid         = fopen(cor_file,'r','native');
%         [rmg,count] = fread(fid,[nx,ny*2],'real*4');
%         status      = fclose(fid);
%         %%mag         = flipud((rmg(1:nx,1:2:ny*2))');
%         cor         = flipud((rmg(1:nx,2:2:ny*2))');
%     end
        cor(cor==0)=NaN;
        data (cor<cor_thresh)= NaN;
        datastruct.data=data;
        datastruct.cor=cor;
end

if(nargin>=9)
    fid  = fopen(demf,'r');
    temp = fread(fid,[nx,ny],'int16');
    datastruct.z=flipud(temp');
    datastruct.hgt=datastruct.z;
end

if(nargin>=12)
    fid  = fopen(conncomp,'r');
    temp = fread(fid,[nx,ny],'uint8');
    conncomp=flipud(temp');
    data(conncomp==0)=NaN;
    if(nargin>=13)
        data(conncomp>conncomp_ind)=NaN;
    end
    %data.conncomp=conncomp;
    datastruct.data=data;
    %%%%    figure;pcolor(conncomp);shading flat;colormap jet;colorbar
end


end
