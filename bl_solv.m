function [int, ils, itr, its, delstar, theta] = bl_solv(x,cp)
global Re duedx ue0
He(1) = 1.57258;    % initialise first He
He_sep = 1.51509;   % store value of He after laminar separation

ue = sqrt(1-cp);

integral = 0;       % initialise momentum integral
n = length(x);      % set maximum index
laminar = true;     % set laminar flag

% initialise placeholders
int = 0;
ils = 0;
itr = 0;
its = 0;

% initialise index, begin calculation for panel
i=1;

% initialise first point
duedx = (ue(i)-0)/(x(i)-0);
integral = integral + ueintbit(0, 0, x(i), ue(i));
theta(i) = sqrt(0.45/Re * ue(i)^(-6) * integral);
m = -Re*(theta(i)^2) * duedx;
H = thwaites_lookup(m);
He(i)=laminar_He(H);
delstar(i) = theta(i)/H;

while laminar && i < n
    i = i + 1;
    duedx = (ue(i)-ue(i-1))/(x(i)-x(i-1));
    integral = integral + ueintbit(x(i-1), ue(i-1), x(i), ue(i));
    theta(i) = sqrt(0.45/Re * ue(i)^(-6) * integral);

    m = -Re*(theta(i)^2) * duedx;
    H = thwaites_lookup(m);
    He(i)=laminar_He(H);
    Rethet = Re.*ue(i).*theta(i);
    delstar(i) = theta(i)/H;

    if log(Rethet) >= 18.4*He(i) - 21.74 %turbulent transition
        laminar = false;
        int=i;
    elseif m >= 0.09  % laminar separation
        laminar = false;
        He(i) = He_sep;
        ils = i;
    elseif laminar ==1 && i == n
        disp(["Did not separate " Re/1e3 x(i) Rethet/1000])
    end   
end

delta(i) = He(end) * theta(end);

while its == 0 && i < n
    i = i+1;
    ue0 = ue(i);
    duedx = (ue(i)-ue(i-1))/(x(i)-x(i-1));
    thick0(1)= theta(i-1);
    thick0(2) = delta(i-1);

    [delx, thickhist]=ode45(@thickdash,[0,x(i)-x(i-1)],thick0);
    
    theta(i) = thickhist(end,1);
    delta(i) = thickhist(end,2);
    
    He(i) = delta(i)/theta(i);
    H=(11*He(i)+15)./(48*He(i)-59);
    delstar(i) = theta(i)/H;
    if  He(i) < 1.46 %turbulent separation
        its = i;
        H=2.803;
    elseif itr == 0 && laminar == false && He(i) > 1.58  % turbulent reattachment
        itr = i;
    end   
end

while its~=0 && i < n
    i = i + 1;
    theta(i) = theta(i-1)*(ue(i-1)/ue(i))^(H+2);
    He(i) = He(i-1);
    delstar(i) = theta(i)/H;
end
end

