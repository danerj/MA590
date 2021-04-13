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
d = 0.1
[Gc,y_true,x_true] = gravity(n,ex,s_a,s_b,d);

K_percent_change = abs(cond(Gc) - cond(G))/cond(G) * 100;

% Modify gravity function to extend length of nonzero source x(t) for t in 
% the interval [-0.5,1.5] (previously, nonzero for t in [0,1])
[~,y_mod,x_mod] = gravity_mod(n,ex,s_a,s_b,d);

[U,S,V] = svd(Gc);
singular_values = diag(S);

%Testing alpha for 100 logarithmically space points between 10^-1 and 10^2.
alphas = logspace(-1,2,100);

x_tikhonov = zeros(n, length(alphas));

for jj = 1:length(alphas)
    
    %The filter factors change for each new choice of alpha. We use equation (1) to find them.
    filter_factors = singular_values.^2 ./ (singular_values.^2 + (alphas(jj))^2);
    
    %The loop below calculates the Tikhonov solution for a given alpha using equation (2). 
    for ii = 1:n
        
        x_tikhonov(:,jj) = x_tikhonov(:,jj) + ...
          filter_factors(ii) * (U(:, ii))' * y_mod * V(:,ii) / singular_values(ii);
      
    end
    
end

error_2norms = zeros(1,length(alphas));

for kk = 1:length(alphas)
    
    error_2norms(kk) = norm(x_tikhonov(:,kk) - x_true);
    
end

%We propose the optimal alpha to be the alpha that minimizes the norm of the residual.
[minnorm, min_index] = min(error_2norms);
good_alpha = alphas(min_index);
minnorm

figure(1)
plot(1:n, x_true, 1:n, x_tikhonov(:,min_index))

figure(2)
plot(error_2norms,alphas)

%b1

y_noise = y_true + 10^-1*randn(length(y_true),1);

for jj = 1:length(alphas)
    
    %The filter factors change for each new choice of alpha. We use equation (1) to find them.
    filter_factors = singular_values.^2 ./ (singular_values.^2 + (alphas(jj))^2);
    
    %The loop below calculates the Tikhonov solution for a given alpha using equation (2). 
    for ii = 1:n
        
        x_tikhonov(:,jj) = x_tikhonov(:,jj) + ...
          filter_factors(ii) * (U(:, ii))' * y_noise * V(:,ii) / singular_values(ii);
      
    end
    
end

error_2norms = zeros(1,length(alphas));

for kk = 1:length(alphas)
    
    error_2norms(kk) = norm(x_tikhonov(:,kk) - x_true);
    
end

%We propose the optimal alpha to be the alpha that minimizes the norm of the residual.
[minnorm, min_index] = min(error_2norms);
good_alpha_b_1 = alphas(min_index);

figure(3)
plot(1:n, x_true, 1:n, x_tikhonov(:,min_index))

%b2

y_mod_noise = y_mod + 10^-1*randn(length(y_true),1);

for jj = 1:length(alphas)
    
    %The filter factors change for each new choice of alpha. We use equation (1) to find them.
    filter_factors = singular_values.^2 ./ (singular_values.^2 + (alphas(jj))^2);
    
    %The loop below calculates the Tikhonov solution for a given alpha using equation (2). 
    for ii = 1:n
        
        x_tikhonov(:,jj) = x_tikhonov(:,jj) + ...
          filter_factors(ii) * (U(:, ii))' * y_mod_noise * V(:,ii) / singular_values(ii);
      
    end
    
end

error_2norms = zeros(1,length(alphas));

for kk = 1:length(alphas)
    
    error_2norms(kk) = norm(x_tikhonov(:,kk) - x_true);
    
end

%We propose the optimal alpha to be the alpha that minimizes the norm of the residual.
[minnorm, min_index] = min(error_2norms);
good_alpha_b2 = alphas(min_index);
minnorm
figure(4)
plot(1:n, x_true, 1:n, x_tikhonov(:,min_index))
