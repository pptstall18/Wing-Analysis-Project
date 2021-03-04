function [lhsmat] = build_lhs_vectorised(xs,ys)
%Builds the matrix of streamfunction contributions from cylinder.
%Matrix [A] in handout, Eqn 6.

np = length(xs)-1;


xa=repmat(xs(1:np),np,1);
ya=repmat(ys(1:np),np,1);
xb=repmat(xs(2:np+1),np,1);
yb=repmat(ys(2:np+1),np,1);
x= xa';
y= ya';

[infa, infb]=panelinf(xa,ya,xb,yb,x,y);

psip(:,1)=infa(:,1);
psip(:,2:np)=infa(:,2:np)+infb(:,1:np-1);
psip(:,np+1)=infb(:,np);

% build lhs matrix
lhsmat = zeros(np+1,np+1);
lhsmat(1,1) = 1;
lhsmat(2:np,:)=psip(2:np,:)-psip(1:np-1,:);
lhsmat(np+1,np+1) = 1;
end

