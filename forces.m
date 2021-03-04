function [cl, cd] = forces(circ, cp, delstarl, thetal, delstaru, thetau)
% Calculates Lift and Drag on an aerofoil

cl = -2*circ;

ue_trailing = sqrt(1-cp(end));
delstar = delstaru(end) + delstarl(end);
theta = thetal(end) + thetau(end);
H_trailing = delstar/theta;

theta_inf = theta * (ue_trailing)^((H_trailing+5)/2);

cd = 2 * theta_inf;
end
