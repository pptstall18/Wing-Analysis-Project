function [xvdv yvdv cpvdv] = vdvfoil ( npan, alpha )
%
%  function [xvdv yvdv cpvdv] = vdvfoil ( npan, alpha )
%
%  Produces coordinates (xvdv,yvdv) and cps for Van de Vooren foil.
%  Geometry parameters: k defines t.e. angle as pi*(2-k)
%                       eps controls thickness
%  NB Panel spacing non-uniform!
%

% parameters for thin (approx 8%) section
k = 1.95;
eps = 0.03;

% % parameters for thick (approx 14%) section
% k = 1.9;
% eps = 0.05;

% panel locations
theta = (0:npan) * 2*pi/npan;
z1 = (cos(theta)-1) + i*sin(theta);
z2 = (cos(theta)-eps) + i*sin(theta);
a = (1+eps)^(k-1) / 2^k;

z = a * z1.^k./z2.^(k-1)  +  1;
xvdv = real(z);
yvdv = imag(z);

% pressure distribution; w is complex potential, zeta the circle plane
dwdzeta = exp(-i*alpha) - exp(i*(alpha - 2*theta)) ...
                                     + 2*i*sin(alpha)*exp(-i*theta);
dzdzeta = (k - (k-1)*z1./z2) .* (z1./z2).^(k-1);
dwdz = dwdzeta(2:npan)./dzdzeta(2:npan);
cpvdv = ones(1,npan+1);
cpvdv(2:npan) = 1 - abs(dwdz).^2;

