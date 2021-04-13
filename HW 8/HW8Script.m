clear all
close all
clc

n = 10^5;

sigma_values = [0.5, 1, 2, 3, 4];
plot_upper_bounds = [5,6,8,15,20];

for k = 1:length(sigma_values)
    
    sigma = sigma_values(k);

    rayleigh_values = rayleigh_generator(sigma,n);
    
    figure(k)
    
    histogram(rayleigh_values,'Normalization','pdf')
    
    x = linspace(0,plot_upper_bounds(k));
    fx = x.*exp(-x.^2 / (2*sigma^2)) / sigma^2;
    
    hold on
    plot(x,fx,'linewidth',2)
    
    xlabel('$x$','interpreter','latex')
    ylabel('$p(x)$','interpreter','latex')
    title(['Rayleigh Probability Density Function, $\sigma = $ ' num2str(sigma)],'interpreter', 'latex')

end