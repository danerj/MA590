function [ X ] = cgls(y,G,k_max)
%
%Inputs: vector y and matrix G of equation Gx = y
%        k_max: the stopping criterion is determined by the user inputting
%        a maximum number of iterations. There is no 'tolerance' criterion.
%Outputs: Matrix X such that the kth column of X contains the k-1 iteration
%of CGLS (the first column corresponds to x_0).

% Initialization

[m,n] = size(G);
X = zeros(n,1);
s = y;
r = G'*y;
p = r;

% The maximum number of iterations for CGLS is n-1, so if the user does not
% specify a maximum iteration, we will just do as many as possible.
if nargin < 3
    k_max = n-1;
end

%CGLS iterations
for k = 0:min(k_max, n-1)
    
    rksquared = (norm(r))^2;
    Gpk = G*p;
    
    alpha = rksquared/(norm(Gpk))^2;
    
    X(:,k+2) = X(:,k+1) + alpha * p;
    s = s - alpha * Gpk;
    r = G'*s;
    
    p = r + (norm(r))^2 / rksquared * p;

end


