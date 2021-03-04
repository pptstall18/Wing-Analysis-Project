%    foil.m
%
%  Script to compare an 2 aerofoil sections
%

addpath('../')
global Re

foil1='foil';
foil2='naca2412';

for i=1:2
    if i==1
        caseref=foil1;
%  Read in the parameter file
parfile = ['Parfiles/' caseref '.txt'];
fprintf(1, '%s\n\n', ['Reading in parameter file: ' parfile])
[section np Re alpha] = par_read(parfile);
section1=section;
Re1=Re;
%  Read in the section geometry
secfile = ['Geometry/' section '.surf'];
[xk yk] = textread ( secfile, '%f%f' );

%  Generate high-resolution surface description via cubic splines
nphr = 5*np;
[xshr yshr] = splinefit ( xk, yk, nphr );

%  Resize section so that it lies between (0,0) and (1,0)
[xsin ysin] = resyze ( xshr, yshr );

%  Interpolate to required number of panels (uniform size)
[xs ys] = make_upanels ( xsin, ysin, np );

%  Assemble the lhs of the equations for the potential flow calculation
A = build_lhswk3 ( xs, ys );
Am1 = inv(A);

%  Loop over alpha values
for nalpha = 1:length(alpha)

%    rhs of equations
  alfrad = pi * alpha(nalpha)/180;
  b = build_rhs ( xs, ys, alfrad );

%    solve for surface vortex sheet strength
  gam = Am1 * b;

%    calculate cp distribution and overall circulation
  [cp circ] = potential_op ( xs, ys, gam );

%    locate stagnation point and calculate stagnation panel length
  [ipstag fracstag] = find_stag(gam);
  dsstag = sqrt((xs(ipstag+1)-xs(ipstag))^2 + (ys(ipstag+1)-ys(ipstag))^2);
  turb_separation_index = 0;
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
  [iunt iuls iutr iuts delstaru thetau] = bl_solv ( su, cpu );

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
  [ilnt ills iltr ilts delstarl thetal] = bl_solv ( sl, cpl );

%    lift and drag coefficients
  [Cl Cd] = forces ( circ, cp, delstarl, thetal, delstaru, thetau );

%    copy Cl and Cd into arrays for alpha sweep plots

  clswp(nalpha) = Cl;
  cdswp(nalpha) = Cd;
  lovdswp(nalpha) = Cl/Cd;

%    screen output

  disp ( sprintf ( '\n%s%5.3f%s', ...
                   'Results for alpha = ', alpha(nalpha), ' degrees' ) )

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

%    save data for this alpha
  fname = ['Data/' caseref '_' num2str(alpha(nalpha)) '.mat'];
  save ( fname, 'Cl', 'Cd', 'xs', 'ys', 'cp','np', ...
         'sl', 'delstarl', 'thetal', 'lowerbl','gam', ...
         'su', 'delstaru', 'thetau', 'upperbl','alfrad',...
         'turb_separation_index', 'Cl_corrected', 'Re')
end

%  save alpha sweep data in summary file

fname = ['Data/' caseref '.mat'];
save ( fname, 'xs', 'ys', 'alpha', 'clswp', 'cdswp', 'lovdswp' )

subplot(2,2,1)
plot(alpha,clswp);
title(['Cl vs \alpha at Re = ',num2str(Re,'%.2g')]);
xlabel("\alpha");
ylabel("Cl");
grid on;
hold on;

subplot(2,2,2)
plot(alpha,cdswp);
title(['Cd vs \alpha at Re = ',num2str(Re,'%.2g')]);
xlabel("\alpha");
ylabel("Cd");
grid on;
hold on;

subplot(2,2,3)
plot(xs,-cp);
title(['Cp vs chord at Re = ',num2str(Re,'%.2g')]);
xlabel('chord');
ylabel('pressure coefficient');
grid on;
hold on;

subplot(2,2,4)
plot(alpha,clswp./cdswp)
title(['^{Cl}/_{Cd} vs \alpha at Re = ',num2str(Re,'%.2g')]);
xlabel("\alpha");
ylabel('^{Cl}/_{Cd}');
grid on;
hold on;
% 
% subplot(2,2,5)
% plot(clswp,cdswp)
% title(['Cd vs Cl at Re = ',num2str(Re,'%.2g')]);
% xlabel("Cl");
% ylabel("Cd");
% grid on;
% hold on;
    else
        caseref=foil2;
%  Read in the parameter file
parfile = ['Parfiles/' caseref '.txt'];
fprintf(1, '%s\n\n', ['Reading in parameter file: ' parfile])
[section np Re alpha] = par_read(parfile);
section2=section;
Re2=Re;
%  Read in the section geometry
secfile = ['Geometry/' section '.surf'];
[xk yk] = textread ( secfile, '%f%f' );

%  Generate high-resolution surface description via cubic splines
nphr = 5*np;
[xshr yshr] = splinefit ( xk, yk, nphr );

%  Resize section so that it lies between (0,0) and (1,0)
[xsin ysin] = resyze ( xshr, yshr );

%  Interpolate to required number of panels (uniform size)
[xs ys] = make_upanels ( xsin, ysin, np );

%  Assemble the lhs of the equations for the potential flow calculation
A = build_lhswk3 ( xs, ys );
Am1 = inv(A);

%  Loop over alpha values
for nalpha = 1:length(alpha)

%    rhs of equations
  alfrad = pi * alpha(nalpha)/180;
  b = build_rhs ( xs, ys, alfrad );

%    solve for surface vortex sheet strength
  gam = Am1 * b;

%    calculate cp distribution and overall circulation
  [cp circ] = potential_op ( xs, ys, gam );

%    locate stagnation point and calculate stagnation panel length
  [ipstag fracstag] = find_stag(gam);
  dsstag = sqrt((xs(ipstag+1)-xs(ipstag))^2 + (ys(ipstag+1)-ys(ipstag))^2);
  turb_separation_index = 0;
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
  [iunt iuls iutr iuts delstaru thetau] = bl_solv ( su, cpu );

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
  [ilnt ills iltr ilts delstarl thetal] = bl_solv ( sl, cpl );

%    lift and drag coefficients
  [Cl Cd] = forces ( circ, cp, delstarl, thetal, delstaru, thetau );

%    copy Cl and Cd into arrays for alpha sweep plots

  clswp(nalpha) = Cl;
  cdswp(nalpha) = Cd;
  lovdswp(nalpha) = Cl/Cd;

%    screen output

  disp ( sprintf ( '\n%s%5.3f%s', ...
                   'Results for alpha = ', alpha(nalpha), ' degrees' ) )

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

%    save data for this alpha
  fname = ['Data/' caseref '_' num2str(alpha(nalpha)) '.mat'];
  save ( fname, 'Cl', 'Cd', 'xs', 'ys', 'cp','np', ...
         'sl', 'delstarl', 'thetal', 'lowerbl','gam', ...
         'su', 'delstaru', 'thetau', 'upperbl','alfrad',...
         'turb_separation_index', 'Cl_corrected', 'Re')
end

%  save alpha sweep data in summary file

fname = ['Data/' caseref '.mat'];
save ( fname, 'xs', 'ys', 'alpha', 'clswp', 'cdswp', 'lovdswp' )

subplot(2,2,1)
plot(alpha,clswp);
set(gca,'Fontn','Times','FontSize',20,'linewidth',1)
legend(section1,section2,'Location','Best');

subplot(2,2,2)
plot(alpha,cdswp);
set(gca,'Fontn','Times','FontSize',20,'linewidth',1)
legend(section1,section2,'Location','Best');

subplot(2,2,3)
plot(xs,-cp);
set(gca,'Fontn','Times','FontSize',20,'linewidth',1)
legend(section1,section2,'Location','Best');

subplot(2,2,4)
plot(alpha,clswp./cdswp);
set(gca,'Fontn','Times','FontSize',20,'linewidth',1)
legend(section1,section2,'Location','Best');


% subplot(2,2,5)
% plot(clswp,cdswp)
% set(gca,'Fontn','Times','FontSize',20,'linewidth',1)
% legend(section1,section2,'Location','Best');
    end
   
end