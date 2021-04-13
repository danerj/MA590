%%

clear all
close all
clc

gauss = @(x) exp(-x.^2/2)/sqrt(2*pi);
alpha = .5;
sigma0 = .08;
sigma1 = .04;
pix = @(x) alpha/sigma0 * gauss(x/sigma0) + (1-alpha)/sigma1 * gauss((x-1)/sigma1);
x = linspace(-0.5,1.5,10^6);
y = pix(x);
figure(1)
plot(x,y,'linewidth', 2)
xlabel('x')
[ymax,ymaxindex] = max(y)
xmax = x(ymaxindex)
hold on
plot(xmax*ones(1,100),linspace(0,ymax),'linewidth', 2)
plot((1-alpha)*ones(1,100), linspace(0,ymax)','linewidth',2)
legend('\pi_{posterior}', 'x_{MAP}', 'x_{CM}','location','northwest')

string1 = '\alpha = '; 'interpreter'; 'latex';
string1 = strcat(string1, string(alpha));
string2 = '\sigma_0 = '; 'interpreter'; 'latex';
string2 = strcat(string2, string(sigma0));
string3 = '\sigma_1 = '; 'interpreter'; 'latex';
string3 = strcat(string3, string(sigma1));
title(['', 'Figure 1', string1, string2, string3, '']);

%\sigma^2 = \alpha \sigma_0^2 + (1-\alpha)(\sigma_1^2 + 1) - (1-\alpha)^2

sigma_squared_1 = alpha*sigma0^2+(1-alpha)*(sigma1^2 + 1)-(1-alpha)^2;

%%
gauss = @(x) exp(-x.^2/2)/sqrt(2*pi);
alpha = .01;
sigma0 = .001;
sigma1 = .1;
pix = @(x) alpha/sigma0 * gauss(x/sigma0) + (1-alpha)/sigma1 * gauss((x-1)/sigma1);
x = linspace(-0.5,1.5,10^6);
y = pix(x);
figure(2)
plot(x,y,'linewidth', 2)
xlabel('x')
[ymax,ymaxindex] = max(y)
xmax = x(ymaxindex)
hold on
plot(xmax*ones(1,100),linspace(0,ymax),'linewidth', 2)
plot((1-alpha)*ones(1,100), linspace(0,ymax)','linewidth',2)
legend('\pi_{posterior}', 'x_{MAP}', 'x_{CM}','location','north')

string1 = '\alpha = '; 'interpreter'; 'latex';
string1 = strcat(string1, string(alpha));
string2 = '\sigma_0 = '; 'interpreter'; 'latex';
string2 = strcat(string2, string(sigma0));
string3 = '\sigma_1 = '; 'interpreter'; 'latex';
string3 = strcat(string3, string(sigma1));
title(['', 'Figure 2', string1, string2, string3, '']);

%\sigma^2 = \alpha \sigma_0^2 + (1-\alpha)(\sigma_1^2 + 1) - (1-\alpha)^2

sigma_squared_2 = alpha*sigma0^2+(1-alpha)*(sigma1^2 + 1)-(1-alpha)^2;

%%
gauss = @(x) exp(-x.^2/2)/sqrt(2*pi);
alpha = .1;
sigma0 = 1;
sigma1 = 1;
pix = @(x) alpha/sigma0 * gauss(x/sigma0) + (1-alpha)/sigma1 * gauss((x-1)/sigma1);
x = linspace(-0.5,1.5,10^6);
y = pix(x);
figure(3)
plot(x,y,'linewidth', 2)
xlabel('x')
[ymax,ymaxindex] = max(y)
xmax = x(ymaxindex)
hold on
plot(xmax*ones(1,100),linspace(0,ymax),'linewidth', 2)
plot((1-alpha)*ones(1,100), linspace(0,ymax)','linewidth',2)
legend('\pi_{posterior}', 'x_{MAP}', 'x_{CM}','location','northwest')

string1 = '\alpha = '; 'interpreter'; 'latex';
string1 = strcat(string1, string(alpha));
string2 = '\sigma_0 = '; 'interpreter'; 'latex';
string2 = strcat(string2, string(sigma0));
string3 = '\sigma_1 = '; 'interpreter'; 'latex';
string3 = strcat(string3, string(sigma1));
title(['', 'Figure 3', string1, string2, string3, '']);


%%
gauss = @(x) exp(-x.^2/2)/sqrt(2*pi);
alpha = .9;
sigma0 = 1;
sigma1 = 1;
pix = @(x) alpha/sigma0 * gauss(x/sigma0) + (1-alpha)/sigma1 * gauss((x-1)/sigma1);
x = linspace(-0.5,1.5,10^6);
y = pix(x);
figure(4)
plot(x,y,'linewidth', 2)
xlabel('x')
[ymax,ymaxindex] = max(y)
xmax = x(ymaxindex)
hold on
plot(xmax*ones(1,100),linspace(0,ymax),'linewidth', 2)
plot((1-alpha)*ones(1,100), linspace(0,ymax)','linewidth',2)
legend('\pi_{posterior}', 'x_{MAP}', 'x_{CM}')

string1 = '\alpha = '; 'interpreter'; 'latex';
string1 = strcat(string1, string(alpha));
string2 = '\sigma_0 = '; 'interpreter'; 'latex';
string2 = strcat(string2, string(sigma0));
string3 = '\sigma_1 = '; 'interpreter'; 'latex';
string3 = strcat(string3, string(sigma1));
title(['', 'Figure 4', string1, string2, string3, '']);
