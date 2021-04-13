clear all
close all
clc

n = 10^5;

sigma_values = [0.5, 1, 2, 3, 4];
plot_upper_bounds = [5,6,8,15,20];

for k = 1:length(sigma_values)
    
    sigma = sigma_values(k);

    figure(k)
    
    x = linspace(0,plot_upper_bounds(k));
    fx = x.*exp(-x.^2 / (2*sigma^2)) / sigma^2;
    
    f_max = sigma*exp(-1/2) / sigma^2;
    
    hold on
    plot(x,fx,'linewidth',2)
    plot(sigma,f_max,'*','linewidth', 2)
    text(sigma,f_max,'($\sigma$, $\pi(\sigma))$', 'interpreter', 'latex')
    plot(sigma*sqrt(pi/2)*ones(100,1),linspace(0,f_max))
    text(sigma*sqrt(pi/2),0,'$E[x] = \sigma \sqrt{\pi/2}$', 'interpreter','latex')
    
    xlabel('$x$','interpreter','latex')
    ylabel('$p(x)$','interpreter','latex')
    title(['Rayleigh Probability Density Function, $\sigma = $ ' num2str(sigma)],'interpreter', 'latex')

end
