
%% Section 1: Using 100 sensors.

clear all
close all
clc

n = 100; %We simulate with 100 seismic sensors

delta_z = 20/n; %The equally spaced sensors go to a depth of 20 meters.

z = (delta_z/2:delta_z:20-delta_z/2)'; %The midpoints between sensors.

G = 0.2*tril(ones(n,n)); %Matrix for midpoint approximation of integral.

s_true = 1./(1000+40*z); %Slowness values based on linear gradient model.

sensor_depths = z+delta_z/2; %Sensors located at z = 0.2,0.4,...,20 meters.

y = log((25+sensor_depths)/25)/40; %Noiseless travel time predictions.

s = G\y; %Slowness value calculations based on synthetic data.

noiseless_error_norm = ...
    norm(s-s_true) %Compare slowness predictions for noiseless data vs. model.
plot1 = plot(z,s,'+', z, s_true,'o')
%make the plot prettier.
xlabel('depth : $$z$$ (meters)','FontSize',16,'interpreter','latex');
ylabel('slowness : $$s,s_{true}$$ (seconds)','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'Figure 1: Slowness Calculations with 100 Sensors (Noiseless Case)' ''},'Interpreter','latex')
names = {'s', 's_{true}'}
legend(plot1,names)

noise = 0.05*10^-3*randn(n,1); %Noise vector with mean 0 and std dev 0.05 milliseconds.
s_noise = G\(y+noise); %Slowness calculation with noisy data.
noisy_error_norm = ...
    norm(s_noise-s_true)%Compare slowness predictions for noisy data vs. model.
figure(2)
plot2 = plot(z,s_noise,'+', z, s_true,'o')
%make the plot prettier.
xlabel('depth : $$z$$ (meters)','FontSize',16,'interpreter','latex');
ylabel('slowness : $$s,s_{true}$$ (seconds)','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'Figure 2: Slowness Calculations with 100 Sensors (Noisy Case)' ''},'Interpreter','latex')
names = {'s', 's_{true}'}
legend(plot2,names)

%See how much adding noise to the data affects solution error.
relative_error_n100 = 100*(norm(s_noise-s_true)-norm(s-s_true))/norm(s-s_true)

%Consider the condition of the matrix used to solve the system in
%understanding how much adding noise may affect results. 
cond(G)

%% Section 2: Using 4 Sensors

clear all
close all
clc

n = 4; %We simulate with 100 seismic sensors

delta_z = 20/n; %The equally spaced sensors go to a depth of 20 meters.

z = (delta_z/2:delta_z:20-delta_z/2)'; %The midpoints between sensors.

G = 5*tril(ones(n,n)); %Matrix for midpoint approximation of integral.

s_true = 1./(1000+40*z); %Slowness values based on linear gradient model.

sensor_depths = z+delta_z/2; %Sensors located at z = 5,10,15...,20 meters.

y = log((25+sensor_depths)/25)/40; %Noiseless travel time predictions.

s = G\y; %Slowness value calculations based on synthetic data.

noiseless_error_norm = ...
    norm(s-s_true) %Compare slowness predictions for noiseless data vs. model.
plot1 = plot(z,s,'+', z, s_true,'o')
%make the plot prettier.
xlabel('depth : $$z$$ (meters)','FontSize',16,'interpreter','latex');
ylabel('slowness : $$s,s_{true}$$ (seconds)','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'Figure 3: Slowness Calculations with 4 Sensors (Noiseless Case)' ''},'Interpreter','latex')
names = {'s', 's_{true}'}
legend(plot1,names)

noise = 0.05*10^-3*randn(n,1); %Noise vector with mean 0 and std dev 0.05 milliseconds.
s_noise = G\(y+noise); %Slowness calculation with noisy data.
noisy_error_norm = ...
    norm(s_noise-s_true)%Compare slowness predictions for noisy data vs. model.
figure(2)
plot2 = plot(z,s_noise,'+', z, s_true,'o')
%make the plot prettier.
xlabel('depth : $$z$$ (meters)','FontSize',16,'interpreter','latex');
ylabel('slowness : $$s,s_{true}$$ (seconds)','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'Figure 4: Slowness Calculations with 4 Sensors (Noisy Case)' ''},'Interpreter','latex')
names = {'s', 's_{true}'}
legend(plot2,names)

%See how much adding noise to the data affects solution error.
relative_error_n4 = 100*(norm(s_noise-s_true)-norm(s-s_true))/norm(s-s_true)

%Consider the condition of the matrix used to solve the system in
%understanding how much adding noise may affect results. 
cond(G)
