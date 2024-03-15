function stats = plotVenusHistogram(nD,Dlower,Dupper,nV,Vlower,Vupper,nAngles,AtmosphereMultiplier,SurfaceAge)

defval('SurfaceAge',700)
stats=[];
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

try
    load(filename,'FLAVRS','Diameters','Velocities','Angles','CraterDiameters','stats');
    VenusCraterHistogram(FLAVRS,Diameters,Velocities,Angles,...
        CraterDiameters,SurfaceAge)
catch
    parallel=1;
    stats=VenusCraterSimulation(nD,Dlower,Dupper,...
        nV,Vlower,Vupper,nAngles,AtmosphereMultiplier,parallel);
end

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