function [lhsmat] = build_lhswk3(xs,ys)
%Builds the matrix of streamfunction contributions from cylinder.
%Matrix [A] in handout, Eqn 6.

np = length(xs)-1;

%initialise psi matrix
psip = zeros(np, np+1);

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

% %compute psi contributions at each point
% for i=1:np
%     for j =1:np+1
%         % special case: first col
%         if j==1
%             [infa1, infb1] = ...
%                 panelinf(xs(j),ys(j),xs(j+1),ys(j+1),xs(i), ys(i));
%             psip(i,j) = infa1;
%         
%         % special case: last col
%         elseif j == np+1
%             [infa1, infb1] = ...
%                 panelinf(xs(j-1),ys(j-1),xs(j),ys(j),xs(i), ys(i));
%             psip(i,j) = infb1;
% 
%         else
%             %evaluate f(a) at j
%             [infa1, infb1] = ...
%                 panelinf(xs(j),ys(j),xs(j+1),ys(j+1),xs(i), ys(i));
%             
%             %evaluate f(b) at j-1
%             [infa2, infb2] = ...
%                 panelinf(xs(j-1),ys(j-1),xs(j),ys(j),xs(i), ys(i));
%             psip(i,j) = infa1 + infb2;
%         end
%     end   
% end

% build lhs matrix
lhsmat = zeros(np+1,np+1);
i = 2:np;
lhsmat(i,:) = psip(i,:)-psip(i-1,:);

% special cases gamma 1 and gamma 0
% modified kutta condition
lhsmat(1,1) = 1;
lhsmat(1,2) = -1;
lhsmat(1,3) = 0.5;
lhsmat(1,np) = 1;
lhsmat(1,np-1) = -0.5;

lhsmat(np+1,np+1) = -1;
lhsmat(np+1,2) = -1;
lhsmat(np+1,3) = 0.5;
lhsmat(np+1,np) = 1;
lhsmat(np+1,np-1) = -0.5;
end


