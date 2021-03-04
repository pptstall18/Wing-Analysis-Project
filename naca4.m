function [xk yk] = naca4naca(id1, id2, id3)
%
%  function [xk yk] = naca4(id)
%
%  Generates 4 digit naca foil according to specification in id.  nc is the
%  number of camber-line points used in generating the surface coordinates
%  (xk,yk).  Points are ordered anti-clockwise from trailing edge, which 
%  is repeated, so length(xk) is 2*(nc-1) + 1.  Note that coordinates are
%  intended to form 'knots' in spline fit; hence fixed, and relatively
%  coarse, discretisation.
%

% geometry parameters from section id designator
m = 0.01 * (id1);
p = 0.1 * (id2);
tmax = 0.01 * (id3);

% camber line and thickness
xc = [0 .0015 .003 0.006 .012 .025:.025:.1 .15:.05:.3 .4:.1:.9 .95 1];
nc = length(xc);
yc(1) = 0;
t(1) = 0;
for ic = 2:nc-1
  if xc(ic)<=p 
    yc(ic) = m * (2*p - xc(ic)) * xc(ic)/p^2;
  else
    yc(ic) = m * (1 - 2*p + 2*p*xc(ic) - xc(ic)^2)/(1 - p)^2;
  end 
  t(ic) = tmax * (2.969*sqrt(xc(ic)) - 1.26*xc(ic) - 3.516*xc(ic)^2 ...
                                        + 2.843*xc(ic)^3 - 1.036*xc(ic)^4);
end
yc(nc) = 0;
t(nc) = 0;

% surface coordinates

xk(1) = 1;
yk(1) = 0;

theta = atan ( (yc(3:nc) - yc(1:nc-2))./(xc(3:nc) - xc(1:nc-2)) );

xk(nc-1:-1:2) = xc(2:nc-1) - 0.5*t(2:nc-1).*sin(theta);
xk(nc+1:2*nc-2) = xc(2:nc-1) + 0.5*t(2:nc-1).*sin(theta);

yk(nc-1:-1:2) = yc(2:nc-1) + 0.5*t(2:nc-1).*cos(theta);
yk(nc+1:2*nc-2) = yc(2:nc-1) - 0.5*t(2:nc-1).*cos(theta);

xk(2*nc-1) = xk(1);
yk(2*nc-1) = yk(1);

% output .surf file in standard format

nk = 2 * (nc-1);
fname = ['Geometry/naca' num2str(id1) num2str(id2) num2str(id3) '.surf'];

fid = fopen(fname,'w');
fprintf ( fid, '%10.6f %10.6f \n', [xk;yk] );
fclose(fid);
