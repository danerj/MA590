%% Plotting for problem 1b First Order Tikhonov

clear all
close all
clc

n = 64;

%The phillips.m function produces Phillips’ test problem, which is a
%discretization of a Fredholm integral equation of the first kind.
[G, y_true, x_true] = phillips(n);

%We generate noisy data y adding a normally distributed perturbation
%to y_true with zero mean and standard deviation 10?4
y = y_true + 10^-4*randn(length(y_true),1);

%We use 50 logarithmically equally spaced points between 10^-4 and 10^0
alphas = logspace(-4,0);

[L,W] = get_l(n,1);
[U,V,X,C,S] = gsvd(G,full(L));
Z = (inv(X))';
lambdas = sqrt(diag(C'*C));
mus = sqrt(diag(S'*S));
generalized_singular_values = lambdas ./ mus;

x_tikhonov = zeros(n, length(alphas));
%We use 50 logarithmically equally spaced points between 10^-4 and 10^0
alphas = logspace(-4,0);

for jj = 1:length(alphas)
    
    filter_factors = zeros(1,n);
    %The loop below calculates the Tikhonov solution for a given alpha using equation (2). 
    for ii = 1:n
        
        generalized_singular_value = lambdas(ii)/mus(ii);
        filter_factors(ii) = generalized_singular_value^2 / ...
            (generalized_singular_value^2 + alphas(jj)^2);
        if isfinite(filter_factors(ii)) == 0
            filter_factors(ii) = 1;
        end
        
        % Since G is square, we have k = 0
        % If a generalized singular value is not finite, set the corresponding filter factor to 1.

        x_tikhonov(:,jj) = x_tikhonov(:,jj) + ...
          filter_factors(ii) * (U(:, ii))' * y * Z(:,ii) / lambdas(ii);
      
    end
    
    plot1 = plot(1:64, x_tikhonov(:,jj), 1:64, x_true, 1:64, filter_factors, 'g*', 'linewidth', 1);
    xlabel('index $$i$$','FontSize',16,'interpreter','latex');
    ylabel('$$x_{\alpha}, x_{true}$$','FontSize',16,'interpreter','latex');
    title({'' 'True vs. Tikhonov Solution for Regularization Parameter $$\alpha$$ = ' num2str(alphas(jj))},'Interpreter','latex')
    names = {'x_{\alpha}' , 'x_{true}', 'f_i'};
    legend(plot1,names, 'location', 'northwest')
    pause(0.4)
    
end

%% Plotting for problem 1b Second Order Tikhonov

clear all
close all
clc

n = 64;

%The phillips.m function produces Phillips’ test problem, which is a
%discretization of a Fredholm integral equation of the first kind.
[G, y_true, x_true] = phillips(n);

%We generate noisy data y adding a normally distributed perturbation
%to y_true with zero mean and standard deviation 10?4
y = y_true + 10^-4*randn(length(y_true),1);

%We use 50 logarithmically equally spaced points between 10^-4 and 10^0
alphas = logspace(-4,0);

[L,W] = get_l(n,2);
[U,V,X,C,S] = gsvd(G,full(L));
Z = (inv(X))';
lambdas = sqrt(diag(C'*C));
mus = sqrt(diag(S'*S));
generalized_singular_values = lambdas ./ mus;

x_tikhonov = zeros(n, length(alphas));
%We use 50 logarithmically equally spaced points between 10^-4 and 10^0
alphas = logspace(-4,0);

for jj = 1:length(alphas)
    
    filter_factors = zeros(1,n);
    %The loop below calculates the Tikhonov solution for a given alpha using equation (2). 
    for ii = 1:n
        
        generalized_singular_value = lambdas(ii)/mus(ii);
        filter_factors(ii) = generalized_singular_value^2 / ...
            (generalized_singular_value^2 + alphas(jj)^2);
        if isfinite(filter_factors(ii)) == 0
            filter_factors(ii) = 1;
        end
        x_tikhonov(:,jj) = x_tikhonov(:,jj) + ...
          filter_factors(ii) * (U(:, ii))' * y * Z(:,ii) / lambdas(ii);
      
    end
    
    plot1 = plot(1:64, x_tikhonov(:,jj), 1:64, x_true, 1:64, filter_factors, 'g*', 'linewidth', 1);
    xlabel('index $$i$$','FontSize',16,'interpreter','latex');
    ylabel('$$x_{\alpha}, x_{true}$$','FontSize',16,'interpreter','latex');
    title({'' 'True vs. Tikhonov Solution for Regularization Parameter $$\alpha$$ = ' num2str(alphas(jj))},'Interpreter','latex')
    names = {'x_{\alpha}' , 'x_{true}', 'f_i'};
    legend(plot1,names, 'location', 'northwest')
    pause(0.4)
    
end