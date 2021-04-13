clear all
close all
clc

n = 64;
% Initializing Phillips Test Problem
[G, y, x] = phillips(n);
[U,S,V] = svd(G);
singular_values = diag(S);
omega_sup = 2 / S(1,1)^2;


%Test the claim from Problem 1a using various omega choices.
test_omegas = [1, 0.5, 0.25, 10^-1, omega_sup - 10^-4,...
    omega_sup-10^-2, 0.05, 10^-2]; 
%One way of comparing landweber vs svd results stored in error_norms.
error_norms = zeros(1,length(test_omegas));

for test = 1:length(test_omegas)

% Matrix to store Landweber iterations. First column corresponds to x_0
% (here we use zero vector) and the i_th column corresponds to x_{i-1}.
X = zeros(64,1);
omega = test_omegas(test);
iterations = 10;

for k = 1:iterations
    
    % Landweber iteration given by equation (1).
    X(:,k+1) = X(:,k) + omega*G'*(y - G*X(:,k));
    
end

figure(test)
plot(1:n, x, 1:n, X(:,end))

%filter factors as defined by equation (2)
filters_10 = 1 - (1- omega*diag(S).^2).^10;



x_svd = 0;
for i = 1:n
    
    x_svd = x_svd + filters_10(i) * ...
        (U(:,i))' * y / singular_values(i) * V(:,i);
    
end

figure(test)
plot1 = plot(1:n, x, '^', 1:n, X(:,end), 'go',...
    1:n, x_svd, 'm*', 'linewidth',1);
title({['Comparison of True, Landweber, and SVD']...
    ['Solutions of Phillips Test Problem $$\omega = $$'...
    num2str(omega)]}, 'interpreter', 'latex');
xlabel('index $$i$$', 'interpreter', 'latex')
ylabel('$$x, x_{svd}^{(10)}, x_{Landweber}^{(10)}$$',...
    'interpreter', 'latex')
names = {'x', 'x_{svd}^{(10)}', 'x_{Landweber}^{(10)}'};
legend(plot1,names, 'location', 'northeast')

% Note that the 11th column of X is the iteration x_10.
error_norms(test) = norm(x_svd - X(:,end));

end
