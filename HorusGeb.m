% Horus Geb - By Jeff Lee & Peter James, Baylor University, 2024

nD = 1000;
Dlower = 500;
Dupper = 64000;
nV = 47;
Vlower = 10;
Vupper = 56;
nAngles = 90;
AtmosphereMultiplier = 1;


%
%% To plot the histogram with BA:
SurfaceAge=100;
stats = plotVenusHistogram_BA(nD,Dlower,Dupper,nV,Vlower,Vupper,nAngles,AtmosphereMultiplier,SurfaceAge);
stats1=stats(stats(:,12)>0 & ~isnan(stats(:,12)),:);

%% To plot the histogram with LW:
SurfaceAge=800;
stats = plotVenusHistogram_LW(nD,Dlower,Dupper,nV,Vlower,Vupper,nAngles,AtmosphereMultiplier,SurfaceAge);
stats1=stats(stats(:,12)>0 & ~isnan(stats(:,12)),:);

%%
% % To get the crater stats:
% parallel = 1;
% stats = VenusCraterSimulation(nD,Dlower,Dupper,...
%     nV,Vlower,Vupper,nAngles,AtmosphereMultiplier,parallel);