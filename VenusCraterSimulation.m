function stats = VenusCraterSimulation(nD,Dlower,Dupper,...
    nV,Vlower,Vupper,nAngles,AtmosphereMultiplier,parallel,SurfaceAge)

% stats = VENUSCRATERSIMULATION(nD,Dlower,Dupper,nV,Vlower,Vupper,ntheta,parallel)
% 
% This code simulates the formation of craters on Venus from projectiles
% that traverse through Venus's atmosphere
%
% INPUTS:
% nD - number of projectile diameters
% Dlower - lowest projectile diameter
% Dupper - largest projectile diameter
% nV - number of projectile velocities
% Vlower - lowest projectile velocity
% Vupper - highest projectile velocity
% nAngles - number of projectile incidence angles
% AtmosphereMultiplier - fraction of the present-day Venus atmosphere
% parallel - logical value indicating if parallel computing will be used
%
% OUTPUT:
% stats is an array with one row per projectile. The columns are:
%
% e.g., stats = VenusCraterSimulation(30,100,64000,11,10,50,30,1);
%       ...should take ~100 seconds to run

% This borrows from code originally developed by: 
% Marissa Dudek (1), R Shane McGary (2)
% (1) University of North Carolina at Chapel Hill, Geological Sciences
% (2) James Madison University, Geology and Environmental Sciences

defval('parallel',1)
defval('SurfaceAge',700)

[FLAVRS,Diameters,Velocities,Angles]=...
    generate_projectile_parameters(nD,Dlower,Dupper,nV,Vlower,Vupper,nAngles);


dt = 0.2;
tic
if parallel==1
    stats = SolveTrajectoryP(FLAVRS,dt,AtmosphereMultiplier);
else
    stats = SolveTrajectory(FLAVRS,dt,AtmosphereMultiplier);
end
toc

CraterDiameters = stats(:,12);
AtmPerMil = round(AtmosphereMultiplier*1000);
if AtmPerMil<1
    AtmStr = '0000';
elseif AtmPerMil<10
    AtmStr = ['000' int2str(AtmPerMil)];
elseif AtmPerMil<100
    AtmStr = ['00' int2str(AtmPerMil)];
elseif AtmPerMil<1000
    AtmStr = ['0' int2str(AtmPerMil)];
else
    AtmStr = int2str(AtmPerMil);
end
filename = ['VenusCraterSimulation_' int2str(nD) '_' int2str(Dlower) '_' ...
    int2str(Dupper) '_' int2str(nV) '_' int2str(Vlower) '_' ...
    int2str(Vupper) '_' int2str(nAngles) '_' AtmStr '.mat'];
save(filename)

%%

VenusCraterHistogram(FLAVRS,Diameters,Velocities,Angles,CraterDiameters,SurfaceAge)


end


function defval(name,value)
if ~ischar(name)
  error(sprintf(['The first argument of DEFVAL ',...
		'has to be a string with a variable name']));
end
si=1;
if evalin('caller',[ 'exist(''' name ''',''var'')'])
  si=evalin('caller',[ 'isempty(' name ')']);
end
if si
  assignin('caller',name,value);
end
end