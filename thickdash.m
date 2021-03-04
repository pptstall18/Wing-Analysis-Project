function [dthickdx] = thickdash(xmx0,thick)
% Calculates the derivative of thick
% where thick is a vector of [theta ; de]

% process global variables
global Re ue0 duedx
ue = ue0 + duedx * xmx0;
Rethet = Re * ue .* thick(1);

He = thick(2)./thick(1); % He = de/theta

% set H according to current He.
if He >= 1.46
    H = (11*He+15)./(48*He-59);
elseif He<1.46
    H=2.803;
end
    
cf = 0.09146*(((H-1).*Rethet).^-0.232).*exp(-1.260*H);
cdiss = 0.010018*((H-1).*Rethet).^(-1/6);

dthickdx = [cf./2 - (H + 2)./ue*duedx .* thick(1);...
    cdiss - 3./ue * duedx .* thick(2)];
end

