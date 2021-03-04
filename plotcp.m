%test for separation - cl curve
clear

%plot aerofoil
plot(xs,ys,"-")

dx = gradient(xs); dy = gradient(ys);
% Normalise both dx and dy by ds
ds = sqrt(dx.^2+dy.^2);
dx = dx./ds ; dy = dy./ds;

dx = abs(cp).*dx ; dy = abs(cp).*dy ;
% We can use a Scaling Factor for the length of Arrow Vector
% Typically, I prefer 0.5
scale = 0.5;

% Determine the point x1,y1 representing the other end of Arrow Vector
x1 = x + scale*dy ;
y1 = y - scale*dx ;

% Find indices of cp less than 0
ind1 = cp<0;
% Set a broader limit for axis, to avoid any warning/error by arrow command
axis([-5 5 -5 5]);
% Plot all negative cp using blue color
arrow([x(ind1) y(ind1)],[x1(ind1),y1(ind1)],'Length',.5,'Width',0.01,'EdgeColor','b','FaceColor','b');
% Plot all positive cp using red color
arrow([x1(~ind1) y1(~ind1)],[x(~ind1),y(~ind1)],'Length',.5,'Width',0.01,'EdgeColor','r','FaceColor','r');
