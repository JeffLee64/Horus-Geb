function stats = SolveTrajectoryP(FLAVRS,dt,AtmosphereMultiplier)
% INPUT:
% FLAVRS - parameter values for each projectile
% dt - the length of a time step (seconds)
% AtmosphereMultiplier - scale the density and pressure of the atmosphere
%                        by a multiplicative constant
%
% OUTPUT:
% stats - an array of projectile and impact parameters

duration = 200;
dim = size(FLAVRS);
nprojectiles = dim(1,1);

ntimesteps = duration/dt;

% Choose any or no atmospheric effects
% ... ['0' off \ '1' on]
% ... ATMOS atmospheric effects
% ... BREAK effect of breakup on the projectile
flag.ATMOS(1) = 1;
flag.BREAK(1) = 1;
% Flag options
% ... vel change in velocity equation 
%     [1958, 1971, 1980, 1995, 2003, 2004, 2005]
% ... shape projectile shape factor [1.21, 1.919]
% ... drag drag coefficient [0.5, 1, 2, 1971, 2016]
% ... ang change in angle equations [1971, 1980, 1993]
% ... mass change in mass equation [1958, 1971, 1980, 1981]
flag.vel(1) = 1971; flag.shape(1) = 1.21; flag.drag(1) = 1971;...
    flag.ang(1) = 1971; flag.mass(1) = 1971;...
    flag.zrho(1) = 1;
    
% Planetary Target Variables

% Call target variables
% ... g gravity m/s2
g.target = 8.87;
% ... H scale height m
H.target = 15900;
% ... Dsc s-c crater m
Dsc.target = 40000; % m, Schaber et al. (1992) - this reference states 40 km for Type D craters (i.e., double ring craters, see p. 13)
% ... rho_t surf. density kg/m3
rho_t = 2830;
% ... diam_lim ab/frag limit m
diam_lim = 1e20; %10000;
max_nfrags = 1e9; % Most number of fragments allowed
SurfaceAirDensity = 67*AtmosphereMultiplier;
vCrack = 200;
KarmanLineAltitude = 200000; % Maximum altitude

% Initial Conditions
% Define global values and matrices to be populated, in alphabetical order
% ... KarmanLineAltitude m maximum altitude
% ... alt m altitude (z)
% ... ang deg. angle of trajectory relative to horizontal
% ... dist m distance (x)
% ... drag - drag coefficient
% ... frag # number of projectile fragments
% ... mass kg projectile mass
% ... massr - ratio of projectile mass to initial mass
% ... rad m radius of projectile fragments
% ... time s travel time
% ... vel m/s velocity
% ... xarea m2 effective cross-sectional area
% ... zrho kg/m3 atmospheric density

% Initialize vectors for all projectiles
Yparent = NaN(nprojectiles, 1); % Has the projectile become lost?
proj_lost = Yparent; % Has the projectile become lost?
proj_impact = Yparent; % Has the projectile impacted?
impact_time_index = Yparent; % Time of impact
lost_time_index = Yparent; % Time when impactor is lost

% minimum velocity needed for crater formation
min_vel = 3000; % meters per second... subsonic?

alpha = 0.1; % Strength scaling parameter.  Acceptable values: 0.05-1.0
CompleteFragmentation = zeros(nprojectiles,1); % This will change to 1 if Ps(j) > YCrit
FragmentationStarted = zeros(nprojectiles,1);
FragmentationTime = NaN(nprojectiles,1);
CrackPropagationTime = NaN(nprojectiles,1);

YCrit = 330e6; % 330 MPa.  This is the maximum strength of a fragment
fragspergen = 2; % Into how many pieces does each fragment subdivide?

% =========================================================
% Crater formation
% Set global values and matrices to be populated, in alphabetical order
% ... (1) column # time step of altitude = 0
% ... (2) ang deg. angle of trajectory relative to horizontal
% ... (3) dist m distance (x)
% ... (4) frag # number of projectile fragments
% ... (5) mass kg projectile mass
% ... (6) massr % projectile mass ratio to initial
% ... (7) vel m/s velocity
% ... (8) rad m projectile fragment radius
% ... (9) xarea m2 effective cross sectional area
% ... (10) time s time to reach surface
impact = zeros(nprojectiles, 11);
% Set variables to be populated in 'crater' matrix
% ... (1) Dt m transient crater diameter
% ... (2) SCid ID simple or complex crater ID (1/simple, 2/complex)
% ... (3) Df m final crater diameter
% ... (4) df m final crater depth
% ... Calcualte using variables from 'impact' matrix
% ...... (2) ang deg. angle of trajectory relative to horizontal
% ...... (4) frag # number of projectile fragments
% ...... (7) rad m projectile fragment radius
% ...... (8) vel m/s velocity
crater = zeros(nprojectiles, 5);


% Atmospheric Effects
% Calculate atmospheric effects
% =========================================================

% f = waitbar(0,'Please wait...');
fprintf('Progress:\n');
% Aim for 100 steps
interval = floor(nprojectiles/100);
fprintf(['\n' repmat('.',1,100) '\n\n']);

% THIS IS THE LOOP THROUGH PROJECTILES & TIME STEPS
parfor k = 1:nprojectiles
    if mod(k,interval)==0
        fprintf('\b|\n'); % Progress bar
    end
    % Calcualte initial conditions
    
    % Initialize vectors for all timesteps
    proj_alt = NaN(1, ntimesteps); proj_ang = proj_alt; 
    proj_dist = proj_alt; proj_frag = proj_alt; proj_fragmass = proj_alt;
    proj_mass = proj_alt; proj_massr = proj_alt; 
    proj_diam = proj_alt; 
    proj_time = proj_alt; 
    proj_vel = proj_alt;
    proj_xarea = proj_alt; proj_xarea_frag = proj_alt; 
    proj_zrho = proj_alt; 
    proj_cloud_radius = proj_alt;
    proj_falls_like_cloud = proj_alt;
    % ... breakup variables
    Ps = proj_alt;
    Yi = proj_alt;
    % populated in order of appearance
    % ---------------------------------------------------------
    % atmosphere
    proj_alt(1) = KarmanLineAltitude;
    % density at altitude (kg/m3)
    % ... see next section for citations and equation validation
    proj_zrho(1) = SurfaceAirDensity * (exp(-proj_alt(1)/H.target));
    % ---------------------------------------------------------
    % projectile characteristics
    % ... angle of trajectroy (deg. from horizontal)
    proj_ang(1) = FLAVRS(k,3);
    % ... distance (x)
    proj_dist(1) = 0;
    % ... number of initial fragments (1)
    proj_frag(1) = 1;
    % ... mass (kg)
    proj_mass(1) = FLAVRS(k,7);
    % ... mass of each fragment (kg)
    proj_fragmass(1) = FLAVRS(k,7);
    % ... mass ratio (%)
    proj_massr(1) = (proj_mass(1)/proj_mass(1))*100;
    % ---------------------------------------------------------
    % projectile area
    proj_xarea(1) = flag.shape(1)*((proj_mass(1)/FLAVRS(k,5))^(2/3));
    % fragment area
    proj_xarea_frag(1) = flag.shape(1)*((proj_mass(1)/FLAVRS(k,5))^(2/3));
    % ---------------------------------------------------------
    % projectile diameter
    proj_diam(1) = (nthroot((FLAVRS(k,7)/((4/3)*(pi)*FLAVRS(k,5))),3))*2;
    proj_cloud_radius(1) = proj_diam(1)/2;
    proj_falls_like_cloud(1)=0;
    % ---------------------------------------------------------
    % projectile velocity
    proj_vel(1) = FLAVRS(k,4);
    % ---------------------------------------------------------
    % Is the projectile lost?
    proj_lost(k) = 0;
    proj_impact(k) = 0;
    
   % ---------------------------------------------------------
    % time
    proj_time(1) = 0;
    % ---------------------------------------------------------
%     Yparent(k) = 10^(2.107 + ( 0.0624*(sqrt(FLAVRS(k,5))) ) );
    switch FLAVRS(k,5)
        case 2700
            Yparent(k) = 1e6;
        case 3300
            Yparent(k) = 1e7;
        case 4500 
            Yparent(k) = 4e7;
        case 7800
            Yparent(k) = 1e8;
        otherwise
            warning('did not recognize projectile density')
    end

for j = 2:ntimesteps
% ---------------------------------------------------------
% Time
proj_time(j) = (j-1)*dt;
% ---------------------------------------------------------
if proj_lost(k)==1 || proj_impact(k)==1
    % Do nothing, except update a few parameters:
%     proj_mass(j) = proj_mass(j-1);
%     proj_alt(j) = proj_alt(j-1);
else % The projectile is NOT lost:

    if proj_alt(j-1) <= 0
        proj_mass(j) = proj_mass(j-1);
        proj_massr(j) = proj_massr(j-1);
    elseif isnan(proj_alt(j-1))
        proj_alt(j-1) = 0;
        proj_alt(j) = 0;
    end 

    % ---------------------------------------------------------
    % Change in atmospheric density (kg/m3)
    if flag.zrho(1) == 1
        proj_zrho(j) = SurfaceAirDensity * (exp(-proj_alt(j-1)/H.target));
    else
        % Herrick and Phillips (pg. 260-261)
        if proj_alt(j-1) >= 70000
            proj_zrho(j) = 0.0789 * (exp(-23.3*((proj_alt(j-1)-70000)/100000)));
        else
            proj_zrho(j) = SurfaceAirDensity * ((1-(proj_alt(j-1)/94800))^5);
        end
    end

    % ---------------------------------------------------------
    % Projectile diameter
    try
        proj_diam(j) = (nthroot((proj_mass(j-1)/((4/3)*(pi)*FLAVRS(k,5))),3))*2;
    catch
        h=warning('something didn''t work')
    end
    if proj_diam(j) <= 0
        proj_diam(j) = 0;
    end

    % ---------------------------------------------------------    
    % FRAGMENTATION
    % ---------------------------------------------------------    
    if FragmentationStarted(k) == 1
        % This is where the cloud disperses
        Vd = proj_vel(j-1)*sqrt(7/2 *(proj_zrho(j)/FLAVRS(k,5)));
        proj_cloud_radius(j) = proj_cloud_radius(j-1) + dt * Vd;
        if isnan(proj_cloud_radius(j))
            warning('this should not happen')
        end
    else
        proj_cloud_radius(j) = proj_cloud_radius(j-1);
    end
    % ---------------------------------------------------------
    % Break up
    if flag.BREAK == 1
        % For Collins, Melosh, and Marcus
        if proj_diam(j) <= diam_lim(1)
            % Set initial frags
            proj_frag(j) = proj_frag(j-1);
            % Yield strength & fragmentation
            % ... Stagnation pressure
            % ... Yield strength of impactor (equation works for 1000-8000 kg/m3)
            Ps(j) = (proj_zrho(j)*(proj_vel(j-1)^2));
            
%             Yi(j) = 10^(2.107 + ( 0.0624*(sqrt(FLAVRS(k,5))) ) );
            % Our strength: dependent on the number of fragments

            Yi(j) = Yparent(k)*(proj_frag(j-1))^alpha;

            % Breakup if stagnation pressure exceeds yield strength
            if Ps(j) > Yi(j)
                if isnan(FragmentationTime(k))
                    CrackPropagationTime(k) = pi*proj_diam(j)/4/vCrack;
                    FragmentationTime(k) = proj_time(j) + CrackPropagationTime(k);
                end
                if proj_time(j)>FragmentationTime(k) && ...
                        FragmentationStarted(k) == 0
                    FragmentationStarted(k) =1;
                end
                % 1) Find the appropriate number of fragments, proj_frag
                %    -> iteratively increase number of fragments through
                %       generations (assume a factor per generation... 2?)
                % 1a) If Ps exceeds 330, complete fragmentation occurs
                if FragmentationStarted(k) == 1
                    if Ps(j)>YCrit && CompleteFragmentation(k) == 0
                        CompleteFragmentation(k) = 1;
                        nfrags=1e9;
                    elseif CompleteFragmentation(k) == 0
                        nfrags = (Ps(j)/Yparent(k))^(1/alpha);
                        nfrags = fragspergen^ceil(log(nfrags)/log(fragspergen)); 
                        if nfrags>max_nfrags
                            nfrags=max_nfrags;
                        end
                    else % Complete Fragmentation occurred in a previous step
                        nfrags = max_nfrags;
                    end
                    proj_frag(j) = nfrags;
                end
                % nfrags is now a power of fragspergen

                % 2) Calculate the altitude of breakup, taking into account
                %    the crack propagation, and time of breakup
                % 2a) Extrapolate back in time to find the precise altitude
                %     and time of breakup
                % 2b) Extrapolate forward in time to account for crack
                %     propagation
                % 3) Calculate velocity of debris cloud expansion
                % 3a) Is it an airburst?  Does that have a different
                %     velocity? Defined as complete fragmentation?
                % 3b) Has the cloud reached its terminal dispersal
                %     velocity?  If so, v_disp = v_disp_terminal
                %       -> this occurs when fragmentation stops
                % 4) Calculate radius of debris cloud

                % 5) Does the projectile behave like a cloud or like
                %    individual fragments?
                %    -> Check to see which has more drag
                %       -> calculating xarea, proj_diameter, effective
                %          density, etc.
                %       -> denote which mode (cloud vs. fragment)

            end % Ps(j) > Yi(j)
        else % proj_diam <= 1km
            % No fragmentation
            proj_frag(j) = proj_frag(j-1);
        end % proj_diam <= 1km
    else
        % if flag.BREAK == 0
        proj_frag(j) = proj_frag(j-1);
    end % flag.BREAK

    % ---------------------------------------------------------
    % CALCULATE THE NEXT TIMESTEP WITH RK4
    % ---------------------------------------------------------    
    if proj_alt(j-1) <= 0
        % If you've impacted the ground, don't worry about updating
        % variables
        proj_impact(k) = 1;
        impact_time_index(k) = j-1;
%         proj_alt(j) = 0;
%         proj_vel(j) = proj_vel(j-1);
%         proj_ang(j) = proj_ang(j-1);
%         proj_mass(j) = proj_mass(j-1);
%         proj_dist(j) = proj_dist(j-1);
    elseif proj_mass(j-1) <=0
        % If projectile is gone, don't worry about updating
        % variables
        proj_lost(k) = 1;
        lost_time_index(k) = (j-1);
%         proj_alt(j) = proj_alt(j-1);
%         proj_vel(j) = proj_vel(j-1);
%         proj_ang(j) = proj_ang(j-1);
%         proj_mass(j) = 0;
%         proj_dist(j) = proj_dist(j-1);
    elseif proj_vel(j-1) < min_vel
        % If the projectile is too slow, we stop updating
        proj_lost(k)=1;
        lost_time_index(k) = (j-1);
    else
        proj_lost(k)=0;
    
        % calculate the acceleration of a particle / acceleration of a
        % cloud
        drag_ratio  = (proj_xarea(j-1)/proj_fragmass(j-1)) / ...
            (pi*proj_cloud_radius(j-1)^2/(proj_mass(j-1)));
        if FragmentationStarted(k) && drag_ratio<1
        % if it falls as a cloud
            proj_falls_like_cloud(j)=1;
            Cloud_xarea = pi*proj_cloud_radius(j-1)^2;
            Cloud_mass = proj_mass(j-1);
            yi = [proj_alt(j-1),proj_vel(j-1),proj_ang(j-1),Cloud_mass,proj_dist(j-1)];
            y = VenusRK4(yi,dt,Cloud_xarea);
            proj_alt(j) = y(1);
            proj_vel(j) = y(2);
            proj_ang(j) = y(3);
            proj_mass(j) = y(4);
            proj_dist(j) = y(5);
        else
        % if it falls as individual rocks:
            proj_falls_like_cloud(j)=0;
            yi = [proj_alt(j-1),proj_vel(j-1),proj_ang(j-1),proj_fragmass(j-1),proj_dist(j-1)];
            y = VenusRK4(yi,dt,proj_xarea_frag(j-1));
            proj_alt(j) = y(1);
            proj_vel(j) = y(2);
            proj_ang(j) = y(3);
            proj_mass(j) = y(4)*proj_frag(j-1);
            proj_dist(j) = y(5);
        end

        if sum(isnan(y))>0
            warning(['RK4 is buggy for projectile: ' int2str(k)])
            proj_lost(k)=1;
            proj_vel(j) = 0;
            proj_mass(j) = 0;
            proj_diam(j) = 0;
        end
        if proj_alt(j) < 0; proj_alt(j) = 0; end
        if proj_mass(j) < 0; proj_mass(j) = 0; end
    end
    % mass of an individual fragment
    proj_fragmass(j) = proj_mass(j) / proj_frag(j);
    % Effective area
    proj_xarea(j) = (flag.shape(1)*(proj_frag(j)^(1/3)))*...
        ((proj_mass(j)/FLAVRS(k,5))^(2/3));
    proj_xarea_frag(j) = proj_xarea(j) / proj_frag(j);

end
end % j = 1:ntimesteps % inner loop

impacttemp = NaN(1,11);

% Crater Formation

% for k = 1:nprojectiles
% Pull final projectile characteristics
% --------------------------------------------
    % Determine time step of impact
    % ... if the projectile alt has stagnated
    try
    if proj_alt(1,ntimesteps-1) > 0 || proj_lost(k)==1
        impacttemp(1) = ntimesteps;
    elseif proj_impact(k) == 1
        if proj_vel(1,find(proj_alt(1,:) == 0,1))> min_vel
            impacttemp(1) = find(proj_alt(1,:) == 0,1);
            % Pull impact characteristics
            impacttemp(2) = proj_ang(1,find(proj_alt(1,:) == 0,1)); 
            impacttemp(3) = proj_dist(1,find(proj_alt(1,:) == 0,1)); 
            impacttemp(4) = proj_frag(1,find(proj_alt(1,:) == 0,1)); 
            impacttemp(5) = proj_mass(1,find(proj_alt(1,:) == 0,1)); 
            impacttemp(6) = proj_massr(1,find(proj_alt(1,:) == 0,1)); 
            impacttemp(7) = NaN; 
            impacttemp(8) = proj_vel(1,find(proj_alt(1,:) == 0,1));
            impacttemp(9) = proj_xarea(1,find(proj_alt(1,:) == 0,1));
            impacttemp(10) = proj_time(1,find(proj_alt(1,:) == 0,1));
            impacttemp(11) = proj_diam(1,find(proj_alt(1,:) == 0,1));
        else
            proj_lost(k)=1;
            impacttemp(1) = ntimesteps;
        end
    else
        % This should only be true if altitudes are negative or NaNs
    end
    catch
        warning('error 2')
    end
    % --------------------------------------------

% "crater" columns: 
% (1) Transient diameter
% (2) Is it simple (1) or complex (2)?
% (3) Final crater diameter
% (4) Final crater depth
% (5) Fragmentation started?  Yes (1) or no (0)?

% Calculate transient crater diameter
% ... (1) Dt m transient crater diameter
% ... (1) Dt km transient crater diameter (converted for depth calcs)
% ....... Calculation walkthrough
% ....... Dt = 1.161 * [ (rho_p/rho_t)^(1/3) * L^0.78 * V^0.44 * g^-0.22 * sin(theta)^(1/3) ]
% ....... crater(1) = 1.161 * ( (FLAVRS(k,5)/rho_t)^(1/3) *
% ((impact(7)*2)^0.78) * (impact(8)^0.44) * (g.target^-0.22) * sin(impact(2))^(1/3) );
% calculate transient diameter


try
cratertemp = NaN(1,5);
try
if proj_lost(k)==0 || (proj_alt(1,ntimesteps-1)) <= 0
    if FragmentationStarted(k)==0
    cratertemp(1) = 1.161 * ( ((FLAVRS(k,5)/rho_t)^(1/3)) * (impacttemp(11)^0.78) * (impacttemp(8)^0.44) * (g.target^-0.22) * sind(impacttemp(2))^(1/3) ) ;
    elseif FragmentationStarted(k)==1
        TransientD = 1.161 * ( ((FLAVRS(k,5)/rho_t)^(1/3)) * (impacttemp(11)^0.78) * (impacttemp(8)^0.44) * (g.target^-0.22) * sind(impacttemp(2))^(1/3) ) ;
        CloudD = 2*proj_cloud_radius(1,impacttemp(1));
        if CloudD<=TransientD
            cratertemp(1) = TransientD;
        else % The cloud is bigger than the transient cavity would be
            cratertemp(1) = 0;
            proj_impact(k) = 0;
        end
    else
        warning('Something Went Wrong')
    end
else % The projectile was lost
    cratertemp(1)=0;
end
catch
    warning('error 4')
    disp(proj_alt(1,ntimesteps-1))
end
cratertemp(1) = cratertemp(1)/1000;
% Determine simple or complex crater
% ... SCid = simple or complex crater identification
% ... (2) SCid ID simple or complex crater ID (1/simple, 2/complex)
% Is it simple or complex?
if cratertemp(1) < (Dsc.target(1)/1000)
    cratertemp(2) = 1;
else
    cratertemp(2) = 2;
end % cratertemp(1) < Dsc.target(1)
% Calculate final crater diameter
% ... Df = Calculate transient crater diameter violin plot
% ... (3) Df km final crater diameter
% ... (4) df km final crater depth
% ... D = 1.17 * ( Dt^1.13 / Dsc.target(1)^0.13 )
if cratertemp(2) == 1
    % Simple crater diameter
    cratertemp(3) = (1.25 * cratertemp(1));
    % Simple crater depth
    a = ( cratertemp(1) / (2*sqrt(2)) );
    b = ( 0.07 * ( (cratertemp(1)^4) / (1.25*(cratertemp(1)))^3) );
    c = ( 2.8 * (0.032 * (1.25*(cratertemp(1)))^3) );
    d = ( ( a + ( 0.04 * ( ((cratertemp(1)^4)) / ((1.25*(cratertemp(1)))^3) ) ) ) / ( a * (1.25*(cratertemp(1)))^2 ) );
    cratertemp(4) = a+b-(c*d);
else % cratertemp(2) == 2
    % Complex crater diameter 
    cratertemp(3) = (1.17 * ( cratertemp(1)^1.13 / (Dsc.target(1)/1000)^0.13 ));
    % Complex crater depth
    cratertemp(4) = 1.04*((cratertemp(3))^0.301);
end % cratertemp(2) == 1
cratertemp(5) = FragmentationStarted(k);

crater(k,:) = cratertemp;
impact(k,:) = impacttemp;

catch
    warning('error 3')
end
end % k = 1:nprojectiles % outer loop

% Pull real numbers only
% impact = real(impact);
% crater = real(crater);

% Model Outputs
% =========================================================
% Pull statistics
% Pull initial and final characteristics
% ... (1) rho_p kg/m3 projectile density
% ... (2) init_diam m initial diameter
% ... (3) init_mass kg initial mass
% ... (4) init_vel km/s initial velocity
% ... (5) init_angle deg initial degree
% ... (6) final_diam m final diameter
% ... (7) final_mass kg final mass
% ... (8) impact time
% ... (9) final_vel km/s final velocity
% ... (10) final_angle deg final degree
% ... (11) trans_diam m transient crater diameter
% ... (12) crater_diam m final crater diameter
% ... (13) crater_depth m final crater depth
% ... (14) Fragmentation started
stats = horzcat(FLAVRS(:,5),FLAVRS(:,2),FLAVRS(:,7),FLAVRS(:,4),FLAVRS(:,3),...
impact(:,11), impact(:,5), impact(:,10), impact(:,8)/1000, impact(:,2),...
(crater(:,1)*1000),(crater(:,3)*1000),(crater(:,4)*1000),crater(:,5));



end