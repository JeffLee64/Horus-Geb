function outvec = TrajectoryDerivativesV(alt,vel,ang,mass,xarea)

htc = 0.02; % proj_htc(k,j)
hoa = 5e6; % proj_hoa(k,j)
dragc = 0.47;
gacc = 8.87; % g.target;
Rtarget = 6052e3; % R.target
surface_rho = 67; % zrho.target
ScaleH = 15900; % H.target

zrho = surface_rho * (exp(-alt/ScaleH));

% Calculate the derivatives:
   altp = -  (vel * (sind(ang)));
   velp = - (1/2)*dragc*zrho*(vel^2)*(xarea/mass)+(gacc*cosd(90-ang));
%    angp = cosd(ang)*(vel/(Rtarget+alt) - gacc/vel);
   angp = abs(-sind(( (gacc/vel) * (1-(cosd(90-ang))^2))));  
   massp = - (1/2)*(htc*zrho *(vel^3)) *(xarea/hoa);    
   distp =  (vel*cosd(ang)) / (1+(alt/Rtarget));

outvec = [altp,velp,angp,massp,distp];

if sum(~isreal(outvec))
    warning('Something isn''t right in TrajectoryDerivatives')
    pause
end

end