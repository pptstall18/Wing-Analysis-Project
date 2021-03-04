%  Script for use in studying resolution requirements of panel method 
%  calculation.  To alter incidence, edit 'alpha' below.  To alter
%  Van de Vooren geometry parameters, see vdvfoil.m.

addpath('../')

clear;
close all;

%  free-stream incidence
alpha = 0;
% alpha = 0;

%  Van de Vooren geometry and pressure distribution
npin = 2000;
[xsin ysin cpex] = vdvfoil( npin, alpha );

% nexttile
figure
plot(xsin,ysin)
axis('equal')
xlabel('x/c')
ylabel('y/c')
title('Van de Vooren aerofoil')

disp('Starting 100 panel calculation ...')
np = 100;
[xs ys] = make_upanels( xsin, ysin, np );

A = build_lhswk3 ( xs, ys );
b = build_rhs ( xs, ys, alpha );

gams = inv(A) * b;
xs1 = xs;
ys1 = ys;
cp1 = 1 - gams.^2;

disp('Starting 200 panel calculation ...')
np = 200;
[xs ys] = make_upanels( xsin, ysin, np );

A = build_lhswk3 ( xs, ys );
b = build_rhs ( xs, ys, alpha );

gams2 = inv(A) * b;
xs2 = xs;
ys2 = ys;
cp2 = 1 - gams2.^2;

disp('Starting 400 panel calculation ...')
np = 400;
[xs ys] = make_upanels( xsin, ysin, np );

A = build_lhswk3 ( xs, ys );
b = build_rhs ( xs, ys, alpha );

gams4 = inv(A) * b;
xs4 = xs;
ys4 = ys;
cp4 = 1 - gams4.^2;

disp('Starting 800 panel calculation ...')
np = 800;
[xs ys] = make_upanels( xsin, ysin, np );

A = build_lhswk3 ( xs, ys );
b = build_rhs ( xs, ys, alpha );

gams8 = inv(A) * b;
xs8 = xs;
ys8 = ys;
cp8 = 1 - gams8.^2;

%plot circulations
delx = gradient(xs1);
dely = gradient(ys1);
del = sqrt(delx.^2 + dely.^2);
circ1 = sum(gams.*del);

delx = gradient(xs2);
dely = gradient(ys2);
del = sqrt(delx.^2 + dely.^2);
circ2 = sum(gams2.*del);

delx = gradient(xs4);
dely = gradient(ys4);
del = sqrt(delx.^2 + dely.^2);
circ4 = sum(gams4.*del);

delx = gradient(xs8);
dely = gradient(ys8);
del = sqrt(delx.^2 + dely.^2);
circ8 = sum(gams8.*del);

% nexttile
figure
box on
plot(xsin,-cpex,xs1,-cp1,'--',xs2,-cp2,'-.',xs4,-cp4,'-+',xs8,-cp8,'-x')
xlabel('x/c')
ylabel('-c_p')
title('Van de Vooren c_ps; varying panel size')
legend('exact','100pans','200pans','400pans','800pans','Location','best','Location','northeast')
set(gca,'Fontn','Times','FontSize',30,'linewidth',1)
set(gcf,'position',[160 280 600 460])
set(gcf,'PaperPositionMode','auto')

% % nexttile
% figure
% plot(xsin,-cpex,xs1,-cp1,'--')
% 
% % nexttile
% figure
% plot(xsin,-cpex,xs2,-cp2,'--')
% 
% % nexttile
% figure
% plot(xsin,-cpex,xs4,-cp4,'--')
% 
% % nexttile
% figure
% plot(xsin,-cpex,xs8,-cp8,'--')
