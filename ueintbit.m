function [f] = ueintbit(xa,ua,xb,ub)

umean=(ua+ub)/2;

udiff=abs(ub-ua);

xdiff=xb-xa;

f=(umean.^5+5/6*umean.^3.*(udiff).^2+1/16*umean.*(udiff)^.4).*xdiff;

end

