%%  
%  Script to analyse an aerofoil section using potential flow calculation
%  and boundary layer solver.

clear

addpath('../')
global Re

error = -0.42;

%  Read in the parameter file
caseref = 'lyndon_flap_optimised';
parfile = ['Parfiles/' caseref '.txt'];
fprintf(1, '%s\n\n', ['Reading in parameter file: ' parfile])
[sectionf, np ,Re ,alphaf] = par_read(parfile);

%  Read in the section geometry
secfile = ['Geometry/' sectionf '.surf'];
[xkf ,ykf] = textread ( secfile, '%f%f' );

%  Generate high-resolution surface description via cubic splines
nphr = 5*np;
[xshrf ,yshrf] = splinefit ( xkf, ykf, nphr );

% Resize section so that it lies between (0,0) and (1,0)
[xsin_temp, ysin_temp] = resyze ( xshrf, yshrf );

% Model 5 degree flaps on approach: twist the secion clockwise by 5 degrees and shift.
twist=(20)*pi/180;

xsinf=0.1.* (xsin_temp*cos(twist)+ysin_temp*sin(twist)) + 1;
ysinf=0.1.* (-xsin_temp*sin(twist)+ysin_temp*cos(twist)) - 0.02;

%  Interpolate to required number of panels (uniform size)
[xsf ,ysf] = make_upanels ( xsinf, ysinf, np );

caseref = 'high_speed_lyndon_p_at_low_speed';
parfile = ['Parfiles/' caseref '.txt'];
fprintf(1, '%s\n\n', ['Reading in parameter file: ' parfile])
[section ,np, Re ,alpha] = par_read(parfile);

%  Read in the section geometry
secfile = ['Geometry/' section '.surf'];
[xk ,yk] = textread ( secfile, '%f%f' );

%  Generate high-resolution surface description via cubic splines
nphr = 5*np;
[xshr ,yshr] = splinefit ( xk, yk, nphr );

%  Resize section so that it lies between (0,0) and (1,0)
[xsin ,ysin] = resyze ( xshr, yshr );

%  Interpolate to required number of panels (uniform size)
[xs, ys] = make_upanels ( xsin, ysin, np );


xa=repmat(xs(1:np),np,1);
ya=repmat(ys(1:np),np,1);
xb=repmat(xs(2:np+1),np,1);
yb=repmat(ys(2:np+1),np,1);

xaf=repmat(xsf(1:np),np,1);
yaf=repmat(ysf(1:np),np,1);
xbf=repmat(xsf(2:np+1),np,1);
ybf=repmat(ysf(2:np+1),np,1);

[infa_largeonsmall, infb_largeonsmall] = panelinf(xa, ya, xb, yb, xaf', yaf');
[infa_smallonlarge, infb_smallonlarge] = panelinf(xaf, yaf, xbf, ybf, xa', ya');

psip(:,1)=infa_largeonsmall(:,1);
psip(:,2:np)=infa_largeonsmall(:,2:np)+infb_largeonsmall(:,1:np-1);
psip(:,np+1)=infb_largeonsmall(:,np);
lhsmat_f = zeros(np+1,np+1);
i = 2:np;
lhsmat_f(i,:) = psip(i,:)-psip(i-1,:);

lhsmat_f(1,1) = 1;
lhsmat_f(1,2) = -1;
lhsmat_f(1,3) = 0.5;
lhsmat_f(1,np) = 1;
lhsmat_f(1,np-1) = -0.5;

lhsmat_f(np+1,np+1) = -1;
lhsmat_f(np+1,2) = -1;
lhsmat_f(np+1,3) = 0.5;
lhsmat_f(np+1,np) = 1;
lhsmat_f(np+1,np-1) = -0.5;

psip(:,1)=infa_smallonlarge(:,1);
psip(:,2:np)=infa_smallonlarge(:,2:np)+infb_smallonlarge(:,1:np-1);
psip(:,np+1)=infb_smallonlarge(:,np);
lhsmat_m = zeros(np+1,np+1);
i = 2:np;
lhsmat_m(i,:) = psip(i,:)-psip(i-1,:);

lhsmat_m(1,1) = 1;
lhsmat_m(1,2) = -1;
lhsmat_m(1,3) = 0.5;
lhsmat_m(1,np) = 1;
lhsmat_m(1,np-1) = -0.5;

lhsmat_m(np+1,np+1) = -1;
lhsmat_m(np+1,2) = -1;
lhsmat_m(np+1,3) = 0.5;
lhsmat_m(np+1,np) = 1;
lhsmat_m(np+1,np-1) = -0.5;

%  Assemble the lhs of the equations for the potential flow calculation
Af = build_lhswk3 ( xsf, ysf );
B = lhsmat_m;

Bf = lhsmat_f;
A = build_lhswk3 ( xs, ys );

%% 
for nalpha = 1:length(alphaf)
%    rhs of equations
  alfrad = pi * alphaf(nalpha)/180;
  bf = build_rhs ( xsf, ysf, alfrad );
  b = build_rhs ( xs, ys, alfrad );
  
%    solve for surface vortex sheet strength
  gam = (Bf-Af*inv(B)*A)\(bf-Af*inv(B)*b);
  gamf = inv(B)*b - inv(B)*A*gam + error;

%    calculate cp distribution and overall circulation
  [cp, circ] = potential_op ( xs, ys, gam );

%    locate stagnation point and calculate stagnation panel length
  [ipstag, fracstag] = find_stag(gam);
  dsstag = sqrt((xs(ipstag+1)-xs(ipstag))^2 + (ys(ipstag+1)-ys(ipstag))^2);
  turb_separation_index = 0;
  lower_laminar_sep_index = 0;
%    upper surface boundary layer calc

%    first assemble pressure distribution along bl
  clear su cpu
  su(1) = fracstag*dsstag;
  cpu(1) = cp(ipstag);
  for is = ipstag-1:-1:1
    iu = ipstag - is + 1;
    su(iu) = su(iu-1) + sqrt((xs(is+1)-xs(is))^2 + (ys(is+1)-ys(is))^2);
    cpu(iu) = cp(is);
  end

%    check for stagnation point at end of stagnation panel
  if fracstag < 1e-6
    su(1) = 0.01*su(2);    % go just downstream of stagnation
    uejds = 0.01 * sqrt(1-cpu(2));
    cpu(1) = 1 - uejds^2;
  end

%    boundary layer solver
  [iunt, iuls, iutr ,iuts ,delstaru, thetau] = bl_solv ( su, cpu );

%    lower surface boundary layer calc

%    first assemble pressure distribution along bl
  clear sl cpl
  sl(1) = (1-fracstag) * dsstag;
  cpl(1) = cp(ipstag+1);
  for is = ipstag+2:np+1
    il = is - ipstag;
    sl(il) = sl(il-1) + sqrt((xs(is-1)-xs(is))^2 + (ys(is-1)-ys(is))^2);
    cpl(il) = cp(is);
  end

%    check for stagnation point at end of stagnation panel
  if fracstag > 0.999999
    sl(1) = 0.01*sl(2);    % go just downstream of stagnation
    uejds = 0.01 * sqrt(1-cpl(2));
    cpl(1) = 1 - uejds^2;
  end

%    boundary layer solver
  [ilnt, ills ,iltr ,ilts, delstarl, thetal] = bl_solv ( sl, cpl );

%    lift and drag coefficients
  [Cl, Cd] = forces ( circ, cp, delstarl, thetal, delstaru, thetau );

%    copy Cl and Cd into arrays for alpha sweep plots

  clswp(nalpha) = Cl;
  cdswp(nalpha) = Cd;
  lovdswp(nalpha) = Cl/Cd;

%    screen output
  disp ( sprintf ( '\n%s%5.3f%s', ...
                   'Results for alpha = ', alphaf(nalpha), ' degrees' ) )

  disp ( sprintf ( '\n%s%5.3f', '  Lift coefficient: ', Cl ) )
  disp ( sprintf ( '%s%7.5f', '  Drag coefficient: ', Cd ) )
  disp ( sprintf ( '%s%5.3f\n', '  Lift-to-drag ratio: ', Cl/Cd ) )

  upperbl = sprintf ( '%s', '  Upper surface boundary layer:' );
  if iunt~=0
    is = ipstag + 1 - iunt;
    upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                        '    Natural transition at x = ', xs(is) );
  end
  if iuls~=0
    is = ipstag + 1 - iuls;
    upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                        '    Laminar separation at x = ', xs(is) );
    if iutr~=0
      is = ipstag + 1 - iutr;
      upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                          '    Turbulent reattachment at x = ', xs(is) );
    end
  end
  if iuts~=0
    is = ipstag + 1 - iuts;
    upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                        '    Turbulent separation at x = ', xs(is) );
    turb_separation_index = is;
  end
  upperbl = sprintf ( '%s\n', upperbl );
  disp(upperbl)

  lowerbl = sprintf ( '%s', '  Lower surface boundary layer:' );
  if ilnt~=0
    is = ipstag + ilnt;
    lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                        '    Natural transition at x = ', xs(is) );
  end
  if ills~=0
    is = ipstag + ills;
    lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                        '    Laminar separation at x = ', xs(is) );
    lower_laminar_sep_index = ills;
    if iltr~=0
      is = ipstag + iltr;
      lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                          '    Turbulent reattachment at x = ', xs(is) );
    end
  end
  if ilts~=0
    is = ipstag + ilts;
    lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                        '    Turbulent separation at x = ', xs(is) );
  end
  lowerbl = sprintf ( '%s\n', lowerbl );
  disp(lowerbl)
  
  Cl_corrected = Cl;
  cl_corrected_swp(nalpha) = Cl_corrected;
  lovdswp_corrected(nalpha) = Cl_corrected/Cd;
  
  if turb_separation_index ~= 0
  %correct circulation value
  circ_corrected = 0;
    for ip = turb_separation_index:length(gam)-1
      %cp_corrected(ip) = 1 - gam(ip)^2;

      dels = sqrt((xs(ip+1)-xs(ip))^2 + (ys(ip+1)-ys(ip))^2);
      circ_corrected = circ_corrected  +  dels * (gam(ip)+gam(ip+1))/2;
    end
  
  [Cl_corrected, ~] = forces ( circ_corrected, cp, delstarl, thetal, delstaru, thetau );
  cl_corrected_swp(nalpha) = Cl_corrected;
  lovdswp_corrected(nalpha) = Cl_corrected/Cd;
  end 

%    calculate cp distribution and overall circulation
  [cpf circf] = potential_op ( xsf, ysf, gamf );

%    locate stagnation point and calculate stagnation panel length
  [ipstagf fracstagf] = find_stag(gamf);
  dsstagf = sqrt((xsf(ipstagf+1)-xsf(ipstagf))^2 + (ysf(ipstagf+1)-ysf(ipstagf))^2);
  turb_separation_index_f = 0;
  lower_laminar_sep_index_f = 0;
%    upper surface boundary layer calc

%    first assemble pressure distribution along bl
  clear su cpu
  su(1) = fracstagf*dsstagf;
  cpu(1) = cpf(ipstagf);
  for is = ipstagf-1:-1:1
    iu = ipstagf - is + 1;
    su(iu) = su(iu-1) + sqrt((xsf(is+1)-xsf(is))^2 + (ysf(is+1)-ysf(is))^2);
    cpu(iu) = cpf(is);
  end

%    check for stagnation point at end of stagnation panel
  if fracstagf < 1e-6
    su(1) = 0.01*su(2);    % go just downstream of stagnation
    uejds = 0.01 * sqrt(1-cpu(2));
    cpu(1) = 1 - uejds^2;
  end

%    boundary layer solver
  [iuntf iulsf iutrf iutsf delstaruf thetauf] = bl_solv ( su, cpu );

%    lower surface boundary layer calc

%    first assemble pressure distribution along bl
  clear sl cpl
  sl(1) = (1-fracstagf) * dsstagf;
  cpl(1) = cpf(ipstagf+1);
  for is = ipstagf+2:np+1
    il = is - ipstagf;
    sl(il) = sl(il-1) + sqrt((xsf(is-1)-xsf(is))^2 + (ysf(is-1)-ysf(is))^2);
    cpl(il) = cpf(is);
  end

%    check for stagnation point at end of stagnation panel
  if fracstagf > 0.999999
    sl(1) = 0.01*sl(2);    % go just downstream of stagnation
    uejds = 0.01 * sqrt(1-cpl(2));
    cpl(1) = 1 - uejds^2;
  end

%    boundary layer solver
  [ilntf illsf iltrf iltsf delstarlf thetalf] = bl_solv ( sl, cpl );

%    lift and drag coefficients
  [Clf Cdf] = forces ( circf, cpf, delstarlf, thetalf, delstaruf, thetauf );

%    copy Cl and Cd into arrays for alpha sweep plots
    Cl_combined = Cl + Clf;
    Cd_combined = Cd + Cdf;
    lovd_combined = Cl_combined/Cd_combined;

  clswp_f(nalpha) = Clf;
  cdswp_f(nalpha) = Cdf;
  lovdswp_f(nalpha) = Clf/Cdf;
  
  clswp_combined(nalpha) = Cl_combined;
  cdswp_combined(nalpha) = Cd_combined;
  lovdswp_combined(nalpha) = lovd_combined;

%    screen output

  disp ( sprintf ( '\n%s%5.3f%s', ...
                   'Results for alpha = ', alphaf(nalpha), ' degrees' ) )

  disp ( sprintf ( '\n%s%5.3f', '  Lift coefficient: ', Clf ) )
  disp ( sprintf ( '%s%7.5f', '  Drag coefficient: ', Cdf ) )
  disp ( sprintf ( '%s%5.3f\n', '  Lift-to-drag ratio: ', Clf/Cdf ) )

  upperbl = sprintf ( '%s', '  Upper surface boundary layer:' );
  if iuntf~=0
    is = ipstagf + 1 - iuntf;
    upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                        '    Natural transition at x = ', xsf(is) );
  end
  if iulsf~=0
    is = ipstagf + 1 - iulsf;
    upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                        '    Laminar separation at x = ', xsf(is) );
    if iutrf~=0
      is = ipstagf + 1 - iutrf;
      upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                          '    Turbulent reattachment at x = ', xsf(is) );
    end
  end
  if iutsf~=0
    is = ipstagf + 1 - iutsf;
    upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                        '    Turbulent separation at x = ', xsf(is) );
    turb_separation_index_f = is;
  end
  upperbl = sprintf ( '%s\n', upperbl );
  disp(upperbl)

  lowerbl = sprintf ( '%s', '  Lower surface boundary layer:' );
  if ilntf~=0
    is = ipstagf + ilntf;
    lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                        '    Natural transition at x = ', xsf(is) );
  end
  if illsf~=0
    is = ipstagf + illsf;
    lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                        '    Laminar separation at x = ', xsf(is) );
    lower_laminar_sep_index_f = illsf;
    if iltrf~=0
      is = ipstagf + iltrf;
      lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                          '    Turbulent reattachment at x = ', xsf(is) );
    end
  end
  if iltsf~=0
    is = ipstagf + iltsf;
    lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                        '    Turbulent separation at x = ', xsf(is) );
  end
  lowerbl = sprintf ( '%s\n', lowerbl );
  disp(lowerbl)
  
  Cl_corrected_f = Clf;
  cl_corrected_swp_combined(nalpha) = Cl_corrected_f + Cl_corrected;
  lovdswp_corrected_f(nalpha) = Cl_corrected_f/Cdf;
  
  if turb_separation_index_f ~= 0
  %correct circulation value
  circ_corrected = 0;
    for ip = turb_separation_index_f:length(gamf)-1
      dels = sqrt((xsf(ip+1)-xsf(ip))^2 + (ysf(ip+1)-ysf(ip))^2);
      circ_corrected = circ_corrected  +  dels * (gamf(ip)+gamf(ip+1))/2;
    end
  
  [Cl_corrected_f, ~] = forces ( circ_corrected, cpf, delstarlf, thetalf, delstaruf, thetauf );
  cl_corrected_swp_combined(nalpha) = Cl_corrected_f + Cl_corrected;
  lovdswp_corrected_f(nalpha) = Cl_corrected_f/Cdf;
  end 

%% plot velocity (pressure) distributions
%Define domain
xmin=-0.5;
xmax=2.0;
ymin=-0.5;
ymax=0.5;

%Define no. of grid lines
nx=61;
ny=51;

for i=1:nx
    for j=1:ny
        xm(i,j)=xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j)=ymin + (j-1)*(ymax-ymin)/(ny-1);
        
        % add freestream streamfunction
        psi(i,j)= ym(i,j)*cos(alfrad) - xm(i,j)*sin(alfrad);
        
        % contribution from wing
        for k = 1:np
            [infa, infb] = panelinf(xs(k),ys(k),xs(k+1),ys(k+1),xm(i,j),ym(i,j));
            psi(i,j)= psi(i,j) + gam(k)*infa + gam(k+1)*infb;        
        end
        
        % contribution from flap
        for k = 1:np
            [infa, infb] = panelinf(xsf(k),ysf(k),xsf(k+1),ysf(k+1),xm(i,j),ym(i,j));
            psi(i,j)= psi(i,j) + gamf(k)*infa + gamf(k+1)*infb;        
        end
    end
end

figure
hold all
box on
MAX_PSI = max(psi,[],'all');
MIN_PSI = min(psi, [], 'all');

[gradx, grady] = gradient(psi, (xmax-xmin)/(nx-1), (ymax-ymin)/(ny-1));
velocity2 = -(gradx.^2 + grady.^2); max_vel = max(velocity2,[],'all'); min_vel = min(velocity2, [], 'all');
c2 = min_vel:0.01:max_vel;
pressure = contourf(xm,ym,velocity2,c2,'edgecolor','none','Displayname','Pressure');
colormap(jet);
c = colorbar;
c.Label.String = 'Pressure: 2/\rho * (P-P_0)';
caxis([min_vel max_vel])
t = get(c,'Limits');
set(c,'Ticks',[t(1),t(2)],...
         'TickLabels',{'Low','High'});

%plot streamlines
c = MIN_PSI:0.05:MAX_PSI;
streamlines = contour(xm,ym,psi,c, 'color','k','DisplayName','Streamlines');

% plot aerofoil and flap
aerofoil = plot(xs,ys,"r-");
aerofoil.Annotation.LegendInformation.IconDisplayStyle = 'off';

filler = fill(xs,ys,'w');
filler.Annotation.LegendInformation.IconDisplayStyle = 'off';

aerofoil = plot(xsf,ysf,"r-");
aerofoil.Annotation.LegendInformation.IconDisplayStyle = 'off';

filler = fill(xsf,ysf,'w');
filler.Annotation.LegendInformation.IconDisplayStyle = 'off';

% % plot boundary layer wing
% delx = gradient(xs); dely = gradient(ys); dels = sqrt(dely.^2 + delx.^2);
% boundary_x = [flip(thetau) thetal] .* dely./dels + xs;
% boundary_y = [flip(thetau) thetal] .* -delx./dels + ys;
% boundarylayer = plot(boundary_x,boundary_y,"b--",'DisplayName','Boundary Layer');
% 
% if turb_separation_index ~= 0
%     % plot turbulent separation
%     turb_sep = plot(xs(turb_separation_index),ys(turb_separation_index),'b*','DisplayName', 'Turbulent Separation');
% end
% 
% % plot boundary layer flap
% delx = gradient(xsf); dely = gradient(ysf); dels = sqrt(dely.^2 + delx.^2);
% boundary_x = [flip(thetauf) thetalf] .* dely./dels + xsf;
% boundary_y = [flip(thetauf) thetalf] .* -delx./dels + ysf;
% boundarylayer = plot(boundary_x,boundary_y,"b--",'DisplayName','Boundary Layer');
% 
% if turb_separation_index_f ~= 0
%     % plot turbulent separation
%     turb_sep = plot(xsf(turb_separation_index_f),ysf(turb_separation_index_f),'b*','DisplayName', 'Turbulent Separation');
% end

title(['Flow Visualisation for Aerofoil at \alpha = ',num2str(alfrad/pi*180)])
xlabel('x/c');
ylabel('y/c');
xlim([xmin xmax])
ylim([ymin ymax])
box on
legend('location','southwest');
set(gca,'Fontn','Times','FontSize',16,'linewidth',1)
set(gcf,'position',[644 400 1020 473])
set(gcf,'PaperPositionMode','auto')
end

%% cl_alpha
figure
hold all
box on
plot(alphaf,clswp_combined)
plot(alphaf,cl_corrected_swp_combined)
title(['^{C_l} vs \alpha at Re (',num2str(Re,'%.2g'),')'])
xlabel("\alpha")
ylabel('{C_l}')
legend('C_L\alpha','C_L\alpha corrected','Location','southwest')
grid on;
set(gca,'Fontn','Times','FontSize',16,'linewidth',1)
set(gcf,'position',[160 280 600 460])
set(gcf,'PaperPositionMode','auto')


%% cl_cd
figure
hold all
box on
plot(alphaf,lovdswp_combined)
title(['^{C_l}/_{C_d} vs \alpha at Re (',num2str(Re,'%.2g'),')'])
xlabel("\alpha")
ylabel('^{C_l}/_{C_d}')
legend('C_L/C_D/\alpha','Location','southwest')
grid on;
set(gca,'Fontn','Times','FontSize',16,'linewidth',1)
set(gcf,'position',[160 280 600 460])
set(gcf,'PaperPositionMode','auto')

%% plot both aerofoils
figure
hold all
plot(xs,ys)
plot(xsf, ysf)
box on
title('Aerofoil Shapes')
xlabel("x/c")
ylabel('y/c')
grid on;
set(gca,'Fontn','Times','FontSize',16,'linewidth',1)
legend('Original High-Speed Aerofoil','Flap','Location','southwest')
ylim([-0.05 0.15])
set(gcf,'position',[945 572 886 413])
set(gcf,'PaperPositionMode','auto')