function [psixy] = psipv(xc,yc,Gamma,x,y)
% where (xc, yc) define coordinates of vortex with strength Gamma
% and the influence of the vortex is measured at (x,y)
% returns a scalar value of psi at the measured point (x,y).

r2=(x-xc)^2 + (y-yc)^2;

psixy= -Gamma/(4*pi)*log(r2);