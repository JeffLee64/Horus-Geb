function [FLAVRS,Diameters,Velocities,Angles]=...
    generate_projectile_parameters(nD,Dlower,Dupper,nV,Vlower,Vupper,nAngles)

% Diameters
Step_factor = (Dupper/Dlower)^(1/(nD-1));
Diameters = NaN(nD,1);
if nD>1
    for i=1:nD
        Diameters(i) = Dupper/Step_factor^(nD-i);
    end
else
    Diameters = Dlower;
end

% Velocities
if nV>1
    dV = (Vupper-Vlower)/(nV-1);
    Velocities = Vlower:dV:Vupper; % km/s
else
    Velocities = Vlower;
end

% Angles
if nAngles>1
    dtheta = 90/nAngles;
    Angles = dtheta/2:dtheta:90-dtheta/2;
else
    Angles = 45;
end

Densities  = [2700; 3300; 4500; 7800];

totallength = length(Angles)*length(Velocities)*length(Densities)*length(Diameters);
projectile=NaN(totallength,4);

count=0;
for m=1:length(Angles)
    for l=1:length(Velocities)
        for j=1:length(Densities)
            for i=1:length(Diameters)
                count=count+1;
                projectile(count,:) = [Diameters(i) Densities(j) Velocities(l) Angles(m)];
            end
        end
    end
end

%Populate Population
nProjectiles = size(projectile,1);
FLAVRS = zeros(nProjectiles, 7);

% The columns of FLAVRS:
% (1) F = effective cross-sectional area (m^2)
% (2) L = projectile diameter (m)
% (3) A = angle of trajectory (deg)
% (4) V = velocity at top of atmosphere (m/s)
% (5) R = projectile density (rho)
% (6) S = projectile strength (Pa)
% (7)     projectile mass (kg)
for k = 1:nProjectiles
    % ---------------------------------------------------------
    % (:,2) L = projectile diameter (m)
    FLAVRS(k,2) = projectile(k,1);
    % ---------------------------------------------------------
    % (:,3) a = angle of trajectory
    FLAVRS(k,3) = projectile(k,4);
    % ---------------------------------------------------------
    % (:,4) V = velocity at top of atmosphere (m/s)
    FLAVRS(k,4) = projectile(k,3)*1000;
    % ---------------------------------------------------------
    % (:,5) R = projectile density (rho)
    FLAVRS(k,5) = projectile(k,2);
    % ---------------------------------------------------------
    % (:,6) S = projectile strength
    FLAVRS(k,6) = 10^(2.107+(0.0624*sqrt(FLAVRS(k,5))));
    % ---------------------------------------------------------
    % (:,7) m = projectile mass
    FLAVRS(k,7) =  ((4/3)*(pi)*((projectile(k,1)/2)^3)*projectile(k,2));
    % ---------------------------------------------------------
    % (:,1) F = effective cross-sectional area
    shape.f = 1.21;
    FLAVRS(k,1) = shape.f*((FLAVRS(k,7)/FLAVRS(k,5))^(2/3));
    % ---------------------------------------------------------
end


end