function [rhsvec] = build_rhs(xs,ys, alpha)
%   Build RHS vector of free-stream contributions to streamfunction

np = length(xs)-1;

%initialise RHS column vector
rhsvec = zeros(np+1,1);

%freestream streamfunction at oncoming incidence alpha
psib = ys*cos(alpha)-xs*sin(alpha);

i = 2:np;
rhsvec(i) = psib(i-1)-psib(i);
end

