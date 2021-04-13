%% (a)
clear all
close all
clc

n=32;
[G,y_true,x_true] = shaw(n);
[U,S,V] = svd(G);
singular_values = diag(S);
y = y_true + 10^-3*randn(length(y_true),1); %simulate noisy data.

a = 20; %Number of alpha values we will use.
alphas = logspace(-5, 1, 20); %alphas log equally spaced.

X = zeros(n,length(alphas));

for jj = 1:length(alphas)
    
    %The filter factors change for each new choice of alpha. We use equation (1) to find them.
    filter_factors = singular_values.^2 ./ (singular_values.^2 + (alphas(jj))^2);
    
    %Loop below calculates zero order Tikhonov solution for a given alpha.
    for ii = 1:n
        
        X(:,jj) = X(:,jj) + ...
          filter_factors(ii) * (U(:, ii))' * y * V(:,ii) / singular_values(ii);
      
    end
    
end


surf(X);
xlabel('$$\alpha_j$$', 'interpreter', 'latex')
ylabel('index $$i$$ of $$x$$', 'interpreter', 'latex')
zlabel('$$x_{\alpha_j,i}$$', 'interpreter', 'latex')
title('Tikhonov Solutions of Shaw Test Problem')

figure(2)
plot(1:32, x_true, 1:32, X(:,8))
figure(3)
plot(1:32, x_true, 1:32, X(:,12))
figure(4)
plot(1:32, x_true, 1:32, X(:,16))

%% (b)

clear all
close all
clc

n=32;
[G,y_true,x_true] = shaw(n);
[U,S,V] = svd(G);
singular_values = diag(S);
y = y_true + 10^-3*randn(length(y_true),1); %simulate noisy data.

alphas = logspace(-5, 1, 10^3); %100 alpha vlaues log equally spaced.

X = zeros(n,length(alphas));

for jj = 1:length(alphas)
    
    %The filter factors change for each new choice of alpha. We use equation (1) to find them.
    filter_factors = singular_values.^2 ./ (singular_values.^2 + (alphas(jj))^2);
    
    %Loop below calculates zero order Tikhonov solution for a given alpha.
    for ii = 1:n
        
        X(:,jj) = X(:,jj) + ...
          filter_factors(ii) * (U(:, ii))' * y * V(:,ii) / singular_values(ii);
      
    end
    
end

error_2norms = zeros(1,length(alphas));

for kk = 1:length(alphas)
    
    error_2norms(kk) = norm(X(:,kk) - x_true);
    
end

%We propose the optimal alpha to be the alpha that minimizes the norm of the residual.
good_alpha_index = find(error_2norms == min(error_2norms))
good_alpha = alphas(good_alpha_index)

x_tikhonov_optimal = X(:,good_alpha_index);

plot1 = plot(1:n, x_tikhonov_optimal, 1:n, x_true, 'linewidth', 0.75);
xlabel('index $$i$$','FontSize',16,'interpreter','latex');
ylabel('$$x_{\alpha}, x_{true}$$','FontSize',16,'interpreter','latex');
title({'' 'True vs. Zeroth Order Tikhonov Solution' 'for Regularization Parameter $$\alpha$$ = ' num2str(good_alpha)},'Interpreter','latex')
names = {'x_{\alpha}' , 'x_{true}'};
legend(plot1,names, 'location', 'northeast')

plot(error_2norms, 'm.', 'linewidth', 0.2)
xlabel('$$\alpha_i$$', 'interpreter', 'latex')
ylabel('Error Norm', 'interpreter', 'latex')
title({'Zero Order Tikhonov Approximation Depending on Regularization Parameter' ''},'Interpreter','latex')


x_alpha_norms = zeros(1,length(alphas));
residual_2norms = zeros(1,length(alphas));

%Finding approximate solutions using truncated svd. It is simpler to to
%calculate x_k for k=1,...,64 even though we will only inspect k=2,..62.
for ii = 1:length(alphas)

   x_alpha_norms(ii) = norm(X(:,ii));
   residual_2norms(ii) = norm(G*X(:,ii) - y);
   
end

loglog(residual_2norms, x_alpha_norms, 'g', 'linewidth', 1);
xlabel('$$||r_k||_2$$','FontSize',16,'interpreter','latex');
ylabel('$$||x_{\alpha}||_2$$','FontSize',16,'interpreter','latex');
title({'' 'Figure 5: Finding optimal regularization parameter with an L-curve' ''},'Interpreter','latex')

