% Plotting for problem 1a

clear all
close all
clc

n = 64;

%The phillips.m function produces Phillips’ test problem, which is a
%discretization of a Fredholm integral equation of the first kind.
[G, y_true, x_true] = phillips(n);
[U,S,V] = svd(G);
singular_values = diag(S);

%We generate noisy data y adding a normally distributed perturbation
%to y_true with zero mean and standard deviation 10?4
y = y_true + 10^-4*randn(length(y_true),1);

%We use 50 logarithmically equally spaced points between 10^-4 and 10^0
alphas = logspace(-4,0);

x_tikhonov = zeros(n, length(alphas));

for jj = 1:length(alphas)
    
    %The filter factors change for each new choice of alpha. We use equation (1) to find them.
    filter_factors = singular_values.^2 ./ (singular_values.^2 + (alphas(jj))^2);
    
    %The loop below calculates the Tikhonov solution for a given alpha using equation (2). 
    for ii = 1:n
        
        x_tikhonov(:,jj) = x_tikhonov(:,jj) + ...
          filter_factors(ii) * (U(:, ii))' * y * V(:,ii) / singular_values(ii);
      
    end
    
    
    plot1 = plot(1:64, x_tikhonov(:,jj), 1:64, x_true, 1:64, filter_factors, 'g*', 'linewidth', 1);
    xlabel('index $$i$$','FontSize',16,'interpreter','latex');
    ylabel('$$x_{\alpha}, x_{true}$$','FontSize',16,'interpreter','latex');
    title({'' 'True vs. Tikhonov Solution for Regularization Parameter $$\alpha$$ = ' num2str(alphas(jj))},'Interpreter','latex')
    names = {'x_{\alpha}' , 'x_{true}', 'f_i'};
    legend(plot1,names, 'location', 'northeast')
    pause(0.6)
    
end
