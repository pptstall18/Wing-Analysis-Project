function [infa, infb] = panelinf(xa,ya,xb,yb,x,y)
% Calculate influence coefficients at position (x,y)
% For a panel with starting coordinates (xa, ya)
% and ending coordinates (xb, yb)

delx =xb-xa;
dely =yb-ya;
del = sqrt(delx.^2 + dely.^2);

X = ((x-xa).*delx + (y-ya).*dely)./del;
Y = ((x-xa).*dely - (y-ya).*delx)./del;

[infa, infb] = refpaninf(del,X,Y);

end
