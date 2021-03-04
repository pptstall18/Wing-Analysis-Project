function [infa, infb] = refpaninf_experimental(del,X,Y)
%   Calculates the influence functions of the end circulation on the
%   streamfunction due to the plate.

%   "trap" to avoid float precision errors/overlow
X(abs(X)<1e-8)=1e-8;
Y(abs(Y)<1e-8)=1e-8;


%   compute I0 and I1 as defined in equations 2 and 3 in handout.
%   these represent components of the influence contribution to the overall
%   streamfunction.
I0 = -1/(4*pi)*((X.*log(X.^2+Y.^2))-((X-del).*log((X-del).^2+Y.^2))-2.*del+2.*Y.*(atan(X./Y)-atan((X-del)./Y)));
I1 = 1/(8*pi)*((X.^2 + Y.^2).*log(X.^2 + Y.^2) - ((X-del).^2 + Y.^2).*log((X-del).^2 + Y.^2)-2*X.*del+del.^2);

%   calculate the influence functions a and b
%   these represent the influence of gamma a and gamma b 
%   (sheet strength at ends) on the streamfunction psi.
infa = (((1-X./del).*I0)-I1./del);
infb = (((X./del).*I0)+I1./del);
