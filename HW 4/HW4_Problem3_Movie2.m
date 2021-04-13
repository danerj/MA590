clear all
close all
clc

n = 100;
[G, y_true, x_true] = wing(n);
[U,S,V] = svd(G);
singular_values = diag(S);

alphas = logspace(-12, 1);
x_tikhonov = zeros(n, length(alphas));

for jj = 1:length(alphas)
    
    %The filter factors change for each new choice of alpha. We use equation (1) to find them.
    filter_factors = singular_values.^2 ./ (singular_values.^2 + (alphas(jj))^2);
    
    %The loop below calculates the Tikhonov solution for a given alpha using equation (2). 
    for ii = 1:n
        
        x_tikhonov(:,jj) = x_tikhonov(:,jj) + ...
          filter_factors(ii) * (U(:, ii))' * y_true * V(:,ii) / singular_values(ii);
      
    end
    
    
    plot1 = plot(1:n, x_tikhonov(:,jj), 1:n, x_true, 'linewidth', 1);
    xlabel('index $$i$$','FontSize',16,'interpreter','latex');
    ylabel('$$x_{\alpha}, x_{true}$$','FontSize',16,'interpreter','latex');
    title({'' 'True vs. Tikhonov Solution for Regularization Parameter $$\alpha$$ = ' num2str(alphas(jj))},'Interpreter','latex')
    names = {'x_{\alpha}' , 'x_{true}'};
    legend(plot1,names, 'location', 'northeast')
    pause(0.6)
    
end