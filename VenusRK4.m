function y = VenusRK4(yi,dt,xarea)

yt = yi;
k1 = TrajectoryDerivativesV(yt(1),yt(2),yt(3),yt(4),xarea);
yt = yi + 0.5*dt*k1;
k2 = TrajectoryDerivativesV(yt(1),yt(2),yt(3),yt(4),xarea);
yt = yi + 0.5*dt*k2;
k3 = TrajectoryDerivativesV(yt(1),yt(2),yt(3),yt(4),xarea);
yt = yi + dt*k3;
k4 = TrajectoryDerivativesV(yt(1),yt(2),yt(3),yt(4),xarea);
dydt = (k1 + 2.0*k2 + 2.0*k3 + k4)/6;
y = yi + dt*dydt;

epsilon1 = 0.001;
epsilon2 = 0.01;
epsilon3 = 0.1;
if norm(dydt./yi)>epsilon1
    % Try again with a series of smaller time steps
    if norm(dydt./yi)>epsilon3
        subdivide = 100;
    elseif norm(dydt./yi)>epsilon2
        subdivide = 10;
    else
        subdivide = 2;
    end
    dt2 = dt/subdivide;
    y = yi;
    for i=1:subdivide
        yt = y;
        k1 = TrajectoryDerivativesV(yt(1),yt(2),yt(3),yt(4),xarea);
        yt = y + 0.5*dt2*k1;
        k2 = TrajectoryDerivativesV(yt(1),yt(2),yt(3),yt(4),xarea);
        yt = y + 0.5*dt2*k2;
        k3 = TrajectoryDerivativesV(yt(1),yt(2),yt(3),yt(4),xarea);
        yt = y + dt2*k3;
        k4 = TrajectoryDerivativesV(yt(1),yt(2),yt(3),yt(4),xarea);
        dydt = (k1 + 2.0*k2 + 2.0*k3 + k4)/6;
        y = y + dt2*dydt;
    end
end


end