function ...
    VenusCraterHistogram_LW(FLAVRS,Diameters,Velocities,Angles,CraterDiameters,SurfaceAge)
% HISTOGRAM CREATOR CODE!!
%
% INPUTS:
% FLAVRS - array of input projectile values (length nprojectiles)
% Diameters - short vector of possible Diameters
% Velocities - short vector of possible Velocities
% Angles - short vector of possible Angles
% CraterDiameters - crater diameter (meters) for each projectile
% SurfaceAge - in millions of years

% Calculate number of impacts per unit time, add to crater diameter histogram bins

tic

% RbE = Relative impact flux on Earth relative to the Moon
% RbV = Relative impact flux on Venus relative to the Moon
% LeFeuvre_Wieczorek_2011_Icarus Chronology.pdf - see Table 4.
RbE = 1.62;
RbV = 1.79;
% rhoast = rhop.ch(1);
% Step_Diameter = 1;
YearMult = 1e6; % per this number of years
LandFrac = 1; % fraction of the planet's surface that is land

% For this script to work, Velocity.csv needs to be formatted with
% velocities in the first column (in ascending order).
% VelTab = xlsread('Velocity.csv');
VelTab = xlsread('VelocityProbability.xlsx');
vVec = VelTab(:,1);
PVec = VelTab(:,2);

TotalProb = sum(PVec); % If probabilities don't add to one, make it so.
PVec = PVec / TotalProb;

nbins = 25;
DiameterHistogramLowerLimit=zeros(nbins,1);

for i=1:nbins
DiameterHistogramLowerLimit(i) = sqrt(2)^(i-3.5) * 1000;
end

NHistogram = zeros(size(DiameterHistogramLowerLimit));

% These are the longer lists of parameters for each projectile
D = FLAVRS(:,2);
theta = FLAVRS(:,3);
velocity = FLAVRS(:,4);
density = FLAVRS(:,5);

nprojectiles = length(D);
% Histogram
errorcount=0;
errorcount1=0;
errorcount2=0;
for q=1:nprojectiles

    DiametersIndex = find(Diameters==D(q));
    if isempty(DiametersIndex)
        warning('Could not find the right Diameter')
    end
    if DiametersIndex==1
        binradius = (Diameters(2)-Diameters(1))/2;
        if binradius>Diameters(1)
            Dlower=Diameters(1);
            if errorcount1==0
            warning('You should consider sampling diameters with higher resolution')
            errorcount1=1;
            end
        else
            Dlower = Diameters(DiametersIndex)-binradius;
        end
        Dupper = Diameters(DiametersIndex)+binradius;
    elseif DiametersIndex==length(Diameters)
        binradius = (Diameters(DiametersIndex)-Diameters(DiametersIndex-1))/2;
        Dlower = Diameters(DiametersIndex)-binradius;
        Dupper = Diameters(DiametersIndex)+binradius;
    else
        Dlower = (Diameters(DiametersIndex)+Diameters(DiametersIndex-1))/2;
        Dupper = (Diameters(DiametersIndex)+Diameters(DiametersIndex+1))/2;
    end
    NDupper = ImpactorRatePerYearVenus(Dupper,density(q),RbE,RbV);
    NDlower = ImpactorRatePerYearVenus(Dlower,density(q),RbE,RbV);
    Nsize(q) = NDlower - NDupper;
    NsizeM(q) = Nsize(q)*YearMult;     

    if theta(q)==Angles(1)
        thetalower = 0;
        thetaupper = (Angles(1)+Angles(2))/2;
    elseif theta(q)==Angles(end)
        thetalower = (Angles(end)+Angles(end-1))/2;
        thetaupper = 90;
    else
        a = find(theta(q)==Angles);
        if isempty(a)
            warning('A projectile angle doesn''t match the vector Angles')
        else
            thetalower = (Angles(a)+Angles(a-1))/2;
            thetaupper = (Angles(a)+Angles(a+1))/2;
        end
    end
    
    if ~isreal(CraterDiameters(q))
        warning('Crater diameter is complex... converting to real')
        CDia = real(CraterDiameters(q));
    else
        CDia = CraterDiameters(q);
    end
    if ~isnan(CDia)
        hbin=find(DiameterHistogramLowerLimit>CDia,1);
        if hbin==1
            if CDia~=0
                 errorcount2=errorcount2+1;
            end
        elseif isempty(hbin)
            % This means the Diameter is larger than our largest bin
            % This is probably a model artifact, so don't add to
            % any bin.
%                     NHistogram(end) = NHistogram(end)+NsizeM(q);
            warning('Crater is larger than the largest bin')
            disp(CDia)
            errorcount=errorcount+1;
        else % This means Diameter is not a NaN, not lower than the lowest bin
            NHistogram(hbin-1) = NHistogram(hbin-1) + ...
                NsizeM(q) * PVelocity(velocity(q)/1000,Velocities,PVec,vVec)...
                * Ptheta(thetalower, thetaupper) * Pdensity(density(q));

            if NsizeM(q) * PVelocity(velocity(q)/1000,Velocities,PVec,vVec)...
                * Ptheta(thetalower, thetaupper) * Pdensity(density(q))<0
                warning('negative Probability')
            end
        end
    end
end

if errorcount>0
disp([int2str(errorcount) ' craters were too large, out of ' int2str(nprojectiles)])
end
if errorcount2>0
disp([int2str(errorcount2) ' craters were below the lowest histogram bin, out of ' int2str(nprojectiles)])
end

CraterDatabase = xlsread('Venus Crater Database.xlsx');
VenusD = CraterDatabase(:,5);
EarthCraterDatabase = xlsread('Earth Crater Database.xlsx');
EarthD = EarthCraterDatabase(:,1);

VHistogram = zeros(size(DiameterHistogramLowerLimit));
for p=1:length(VenusD)
    if ~isreal(VenusD(p,1))
        warning('Crater diameter is complex... converting to real')
        VDia = real(VenusD(p,1));
    else
        VDia = VenusD(p,1);
    end
    VDia = VDia*1000;
    if ~isnan(VDia)
%             if Dia == 0
            hbin=find(DiameterHistogramLowerLimit>VDia,1);
            if hbin==1
%                     warning('Diameter is less than lowest histogram limit??')
            elseif isempty(hbin)
                % This means the Diameter is larger than our largest bin
                VHistogram(end) = VHistogram(end)+1;
            else % This means Diameter is not a NaN, not lower than the lowest bin
                VHistogram(hbin-1) = VHistogram(hbin-1) + 1;
            end
    end
end

EHistogram = zeros(size(DiameterHistogramLowerLimit));
for p=1:length(EarthD)
    if ~isreal(EarthD(p,1))
        warning('Crater diameter is complex... converting to real')
        EDia = real(EarthD(p,1));
    else
        EDia = EarthD(p,1);
    end
    EDia = EDia*1000;
    if ~isnan(EDia)
%             if Dia == 0
            hbin=find(DiameterHistogramLowerLimit>EDia,1);
            if hbin==1
%                     warning('Diameter is less than lowest histogram limit??')
            elseif isempty(hbin)
                % This means the Diameter is larger than our largest bin
                EHistogram(end) = EHistogram(end)+1;
            else % This means Diameter is not a NaN, not lower than the lowest bin
                EHistogram(hbin-1) = EHistogram(hbin-1) + 1;
            end
    end
end


figure
% standard deviation for Poisson distribution (times 2 for 95% confidence)
errplus = 2*sqrt(NHistogram*SurfaceAge);
errminus = 2*errplus; 
for i=1:length(errminus)
    if errminus(i)>NHistogram(i)*SurfaceAge
        errminus(i)=NHistogram(i)*SurfaceAge-1e-10;
    end
end
errorbar(DiameterHistogramLowerLimit,NHistogram*SurfaceAge,errminus,errplus,'k')
% bar(DiameterHistogramLowerLimit,NHistogram*SurfaceAge)
set(gca,'YScale','log','XScale','log')
% set(gca,'XScale','log')
hold on

SS = sprintf('Number of Venus Impacts per Bin for a %0.1e Million-', SurfaceAge);
SS = regexprep(SS, 'e\+?(-?\d+)', ' x 10^{$1}');
title([SS, 'Year-Old Atmosphere'])
% xlim([0 20])
% ylim([0 100])
xlabel('Crater Diameter (m)')
ylabel('Number of Impacts')
xlim([5e2 900e3])
% ylim([0 2e2])
ylim([1e-1 1e3])
% somenames={'1-1.1';'1.1-1.2';'1.2-1.3';'1.3-1.4';'1.4-1.5';'1.5-1.6';'1.6-1.7';'1.7-1.8';'1.8-1.9';'1.9-2';'2-3';'3-4';'4-5';'>5'};
% set(gca,'xticklabel',somenames)
hold on
plot([1 1e7],[1 1],':')
plot(DiameterHistogramLowerLimit,VHistogram,'ro-')
plot(DiameterHistogramLowerLimit,EHistogram,'go-')

legend('Expectation','','Venus craters','Earth craters')

% figure
% bar(DiameterHistogramLowerLimit,VHistogram)
% set(gca,'YScale','log','XScale','log')
% SS = sprintf('Number of Land Impacts per Bin per %0.1e', YearMult);
% SS = regexprep(SS, 'e\+?(-?\d+)', ' x 10^{$1}');
% title([SS, ' Years'])
% % xlim([0 20])
% % ylim([0 100])
% xlabel('Crater Diameter')
% ylabel('Number of Impacts')
% xlim([1e3 300e3])



toc

% function ND = ImpactorRatePerYearVenus(D,rhoast,RbE,RbV)
% if (pi*rhoast/6*D^3) > 4e12
%     ND = 10^4.5*(RbE/RbV) * (pi*rhoast/6*D^3)^(-0.75565);
% else
%     ND = (RbE/RbV) * 1e-5;
% end
% end

% function ND = ImpactorRatePerYearVenus(D,~,RbE,RbV)
% ND = (RbV/RbE) * 10^(1.0 + 0.31656*(log10(D/1000)) + 0.10393*(log10(D/1000))^2 + 5.7091e-2*(log10(D/1000))^3 - 8.1475e-2*(log10(D/1000))^4 ...
% - 2.9864e-2*(log10(D/1000))^5 + 1.3977e-2*(log10(D/1000))^6 + 5.8676e-3*(log10(D/1000))^7 - 4.6476e-4*(log10(D/1000))^8 ...
% - 3.8428e-4*(log10(D/1000))^9 - 3.7825e-5*(log10(D/1000))^10);
% end

% function ND = ImpactorRatePerYearVenus(D,rhoast,RbE,RbV) % Bland and Artemieva
% ND = 10^4.5*(RbV/RbE) * (pi*rhoast/6*D^3)^(-0.75565);
% end


% function ND = ImpactorRatePerYearVenus(D,~,RbE,RbV) % LeFeuvre & Wieczorek (1e-4 - 100 km)
% ND = (RbV/RbE) * 1e-6*(D/1000)^(-2.404);
% end


% function ND = ImpactorRatePerYearVenus(D,~,RbE,RbV) % LeFeuvre & Wieczorek (1 - 100 km)
% ND = (RbV/RbE) * 9e-7*(D/1000)^(-1.878);
% end


% function ND = ImpactorRatePerYearVenus(D,~,RbE,RbV) % LeFeuvre & Wieczorek (30 m - 12 km)
% ND = (RbV/RbE) * 9e-7*(D/1000)^(-2.118);
% end


% function ND = ImpactorRatePerYearVenus(D,~,RbE,RbV) % Bland & Artemieva (1 - 100 km)
% ND = (RbV/RbE) * 2e-5*(D/1000)^(-2.267);
% end


% function ND = ImpactorRatePerYearVenus(D,~,RbE,RbV) % Lee (100 m - 64 km)
% ND = (RbV/RbE) * (-0.2693)*(D/1000)^6 + 0.449*(D/1000)^5 + 0.7894*(D/1000)^4 ...
%     - 0.9466*(D/1000)^3 - 0.6397*(D/1000)^2 - 1.6626*(D/1000) - 5.9846;
% end


% function ND = ImpactorRatePerYearVenus(D,~,RbE,RbV) % LeFeuvre & Wieczorek modified (1e-4 - 100 km)
% ND = (RbV/RbE) * 7e-6*(D/1000)^(-1.95);
% end


function ND = ImpactorRatePerYearVenus(D,~,RbE,RbV) % LeFeuvre & Wieczorek (100 m - 64 km)
ND = (RbV/RbE) * 1e-6*(D/1000)^(-1.939);
end


function PT = Ptheta(theta1,theta2)
PT = (cosd(theta1))^2 - (cosd(theta2))^2;
end


function PV = PVelocity(projectileV,Vparameters,tableP,tableV)
    index = find(Vparameters==projectileV);
    if isempty(index)
        warning('velocity is not in the parameter list')
        disp('stop here')
    end
    if index==1
        upperlimit = (Vparameters(1)+Vparameters(2))/2;
        inrange = tableV<upperlimit;
    elseif index==length(Vparameters)
        lowerlimit = (Vparameters(end)+Vparameters(end-1))/2;
        inrange = tableV>lowerlimit;
    else
        upperlimit = (Vparameters(index)+Vparameters(index+1))/2;
        lowerlimit = (Vparameters(index)+Vparameters(index-1))/2;
        inrange = tableV>lowerlimit & tableV<upperlimit;
    end
    PV = sum(tableP.*inrange);
end

function P = Pdensity(rho)
rhop.ch(1) = 2700; %3110, 2700
rhop.ir(1) = 7800; %7600, 7800
rhop.si(1) = 4500; %4396, 4500
rhop.st(1) = 3300; %3320, 3300
    if rho==rhop.ch(1)
        P=0.753; % 100 - 5.7 - 2 - 17 = 75.3%
    elseif rho==rhop.ir(1) 
        P=0.057; % https://en.wikipedia.org/wiki/Iron_meteorite#:~:text=Although%20they%20are%20fairly%20rare,as%20opposed%20to%20stony%20meteorites.
    elseif rho == rhop.si(1)
        P=0.02; % https://geology.com/meteorites/stony-iron-meteorites.shtml#:~:text=Compared%20to%20the%20other%20two,2%25%20of%20all%20known%20meteorites.
    elseif rho == rhop.st(1)
        P=0.17; % https://en.wikipedia.org/wiki/S-type_asteroid#:~:text=S%2Dtype%20asteroids%20are%20asteroids,after%20the%20carbonaceous%20C%2Dtype.
    else
        warning('The projectile''s density was not in our lookup table')
    end
end

end