clear all
close all
clc

% 'gravity' test problem: one-dimensional gravity surveying model problem
% (see Regularization Tools for more details)
n   = 64; % size of problem 
ex  = 3;  % type of source function, here step function
s_a = 0;  % interval of s-axis [s_a,s_b]
s_b = 1;
d   = 0.05; % depth

[G,y_true,x_true] = gravity(n,ex,s_a,s_b,d);

% Modify gravity function to extend length of nonzero source x(t) for t in 
% the interval [-0.5,1.5] (previously, nonzero for t in [0,1])
[~,y_mod,x_mod] = gravity_mod(n,ex,s_a,s_b,d);


[U,S,V] = svd(G);
singular_values = diag(S);

%Testing alpha for 100 logarithmically space points between 10^-1 and 10^2.
%Theses values used because the singular values of G all lie between 10^-1
%and 10^2.
alphas = logspace(-1,2,100);

x_tikhonov = zeros(n, length(alphas));

%Part b part 1

y_noise = y_true + 10^-1*randn(length(y_true),1);

for jj = 1:length(alphas)
    
    %The filter factors change for each new choice of alpha.
    filter_factors = singular_values.^2 ./...
        (singular_values.^2 + (alphas(jj))^2);
    
    %The loop below calculates the Tikhonov solution for a given alpha. 
    for ii = 1:n
        
        x_tikhonov(:,jj) = x_tikhonov(:,jj) + filter_factors(ii) * ...
            (U(:, ii))' * y_noise * V(:,ii) / singular_values(ii);
      
    end
    
end

error_2norms = zeros(1,length(alphas));

for kk = 1:length(alphas)
    
    error_2norms(kk) = norm(x_tikhonov(:,kk) - x_true);
    
end

%We propose the optimal alpha to be the alpha that minimizes the norm of
%the residual.
[minnorm, min_index] = min(error_2norms);
good_alpha = alphas(min_index);

figure(1)
plot1 = plot(1:n, x_true, 1:n, x_tikhonov(:,min_index), 1:n, filter_factors, 'g*', 'linewidth', 1);
xlabel('index $$i$$','FontSize',16,'interpreter','latex');
ylabel('$$x_{\alpha}, x_{true}$$','FontSize',16,'interpreter','latex');
title({'Adding Noise and Committing Crimes:' 'True vs. Zeroth Order Tikhonov Solution' 'for Regularization Parameter $$\alpha$$ = ' num2str(good_alpha)},'Interpreter','latex')
names = {'x_{true}','x_{\alpha}' , 'f_i'};
legend(plot1,names, 'location', 'northeast')

%Calculate the naive solution
x_naive = G \ y_noise;

figure(2)
plot1 = plot(1:n, x_true, 1:n, x_naive,'linewidth', 1);
xlabel('index $$i$$','FontSize',16,'interpreter','latex');
ylabel('$$x_{\alpha}, x_{true}$$','FontSize',16,'interpreter','latex');
title({'Adding Noise and Committing Crimes:' 'True vs. Naive Solution' ''},...
    'Interpreter','latex')
names = {'x_{true}','x_{naive}'};
legend(plot1,names, 'location', 'northeast')

%Part b part 2

y_mod_noise = y_mod + 10^-1*randn(length(y_true),1);
x_tikhonov = zeros(n, length(alphas));

for jj = 1:length(alphas)
    
    %The filter factors change for each new choice of alpha.
    filter_factors = singular_values.^2 ./...
        (singular_values.^2 + (alphas(jj))^2);
    
    %The loop below calculates the Tikhonov solution for a given alpha. 
    for ii = 1:n
        
        x_tikhonov(:,jj) = x_tikhonov(:,jj) + filter_factors(ii) *...
            (U(:, ii))' * y_mod_noise * V(:,ii) / singular_values(ii);
      
    end
    
end

error_2norms = zeros(1,length(alphas));

for kk = 1:length(alphas)
    
    error_2norms(kk) = norm(x_tikhonov(:,kk) - x_true);
    
end

%We propose the optimal alpha to be the alpha that minimizes the norm of
%the residual.
[minnorm, min_index] = min(error_2norms);
good_alpha = alphas(min_index);
figure(3)
plot1 = plot(1:n, x_true, 1:n, x_tikhonov(:,min_index), 1:n, filter_factors, 'g*', 'linewidth', 1);
xlabel('index $$i$$','FontSize',16,'interpreter','latex');
ylabel('$$x_{\alpha}, x_{true}$$','FontSize',16,'interpreter','latex');
title({'Adding Noise and Abiding the Law:' 'True vs. Zeroth Order Tikhonov Solution' 'for Regularization Parameter $$\alpha$$ = ' num2str(good_alpha)},'Interpreter','latex')
names = {'x_{true}','x_{\alpha}' , 'f_i'};
legend(plot1,names, 'location', 'northeast')

%Calculate the naive solution
x_naive = G \ y_mod_noise;

figure(4)
plot1 = plot(1:n, x_true, 1:n, x_naive,'linewidth', 1);
xlabel('index $$i$$','FontSize',16,'interpreter','latex');
ylabel('$$x_{\alpha}, x_{true}$$','FontSize',16,'interpreter','latex');
title({'Adding Noise and Abiding the Law:' 'True vs. Naive Solution' ''},...
    'Interpreter','latex')
names = {'x_{true}','x_{naive}'};
legend(plot1,names, 'location', 'northeast')