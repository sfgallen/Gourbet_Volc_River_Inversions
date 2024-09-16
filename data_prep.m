clear; clc

addpath(genpath('C:\Users\sfgallen\Documents\topo_toolbox\topotoolbox-master'));

% load the prep and data
DEM = GRIDobj('maurice-proj-filled-clip.tif');
pDEM = GRIDobj('pre_inc_maurice16.tif');

DEM.Z(DEM.Z <= 0) = nan;
pDEM.Z(pDEM.Z <= 0) = nan;

DEM.Z(1:23,1:665) = nan; % some odd sliver of data in NW corner removed
pDEM.Z(isnan(DEM.Z)) = nan;

[DEM,pDEM] = largest_overlapping_extent(DEM,pDEM);

DEM = fillsinks(DEM);
pDEM = fillsinks(pDEM);

FD = FLOWobj(DEM);
S = STREAMobj(FD,'minarea',1e6./DEM.cellsize^2);

% picking a location with some incision in the basin and cropping the DEM
xc = -87780;
yc=7738000;
%[x,y] = coord2sub(DEM,xc,yc);
[~,~,ix] = snap2stream(S,xc,yc);

db = drainagebasins(FD,ix);

DEM.Z(db.Z <= 0) = nan;
pDEM.Z(db.Z <= 0) = nan;

DEM = crop(DEM);
pDEM = crop(pDEM);
[DEM,pDEM] = largest_overlapping_extent(DEM,pDEM);

FD = FLOWobj(DEM);
A = flowacc(FD).*DEM.cellsize^2;

% select largest river system for this example
S = STREAMobj(FD,'minarea',1e6./DEM.cellsize^2);
%S = klargestconncomps(S,1);

figure
imageschs(DEM,pDEM); hold on
plot(S,'w-')

Sz = DEM.Z(S.IXgrid);
pSz = pDEM.Z(S.IXgrid);

figure
subplot(2,1,1)
plot(S.distance./1000,Sz,'k.'); hold on
plot(S.distance./1000,pSz,'b.'); hold on
xlabel('Distance (km)'); ylabel('Elevation (m)');
subplot(2,1,2)
plot(S.distance./1000,pSz-Sz,'g.'); hold on
xlabel('Distance (km)'); ylabel('Total Incision (m)');

% make a data structure with relevant data and save it
stream_data.S = S;
stream_data.Sz = Sz;
stream_data.pSz = pSz;
stream_data.S_DA = A.Z(S.IXgrid);

save('stream_data.mat','stream_data');
