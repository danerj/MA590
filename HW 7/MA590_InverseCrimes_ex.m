% MA 590: ST: Computational Inverse Problems (Spring 2020)
% Instructor: Andrea Arnold
%
% Inverse Crimes example script

clear; close all; clc

% 'gravity' test problem: one-dimensional gravity surveying model problem
% (see Regularization Tools for more details)
n   = 64; % size of problem 
ex  = 3;  % type of source function, here step function
s_a = 0;  % interval of s-axis [s_a,s_b]
s_b = 1;
d   = 0.05; % depth

[G,y_true,x_true] = gravity(n,ex,s_a,s_b,d);

cond(G) % check condition of G

% Compute t, s values for solution, data intervals: x = x(t), y = y(s)
dt = (1-0)/n;
ds = (s_b-s_a)/n;
t_vals = 0 + dt*((1:n)' - 0.5);  
s_vals = s_a + ds*((1:n)' - 0.5);

figure(1);
subplot(2,1,1); plot(t_vals,x_true,'-k','LineWidth',2); title('x_{true}'); set(gca,'FontSize',20);
subplot(2,1,2); plot(s_vals,y_true,'-k','LineWidth',2); title('y_{true}'); set(gca,'FontSize',20);
hold on;

% Add noise to data, if desired
noiselev = 1;
y_data = y_true + noiselev*randn(size(y_true));

figure(1);
subplot(2,1,2); plot(s_vals,y_data,'.b','MarkerSize',20); title('y_{true}'); set(gca,'FontSize',20);


% Compute ``naive'' solution via backslash command (proxy for inv(G)*y)
x_sol = G\y_data;

figure(2);
plot(t_vals,x_true,'-k','LineWidth',2); hold on;
plot(t_vals,x_sol,'--r','LineWidth',2); hold off
title('Naive solution');
legend('x_{true}','x_{sol}');
set(gca,'FontSize',20);


%%%%%

% Modify gravity function to extend length of nonzero source x(t) for t in 
% the interval [-0.5,1.5] (previously, nonzero for t in [0,1])
[~,y_mod,x_mod] = gravity_mod(n,ex,s_a,s_b,d);

% Recompute values for t interval: x = x(t)
dt = (1.5+0.5)/n; 
t_mod = -0.5 + dt*((1:n)' - 0.5); 

figure(3);
subplot(2,1,1); plot(t_mod,x_mod,'-k','LineWidth',2); title('x_{true}'); set(gca,'FontSize',20);
subplot(2,1,2); plot(s_vals,y_mod,'-k','LineWidth',2); title('y_{true}'); set(gca,'FontSize',20);
hold on;

% Add noise to data, if desired
y_datamod = y_mod + noiselev*randn(size(y_mod));

figure(3);
subplot(2,1,2); plot(s_vals,y_datamod,'.m','MarkerSize',20); title('y_{true}'); set(gca,'FontSize',20);
hold on;


% Compare solutions, data for previous model to the modified model (with 
% longer nonzero source function)
figure(4);
subplot(2,1,1); plot(t_vals,x_true,'-k','LineWidth',2); title('x_{true}'); set(gca,'FontSize',20);
hold on;
plot(t_mod,x_mod,'--m','LineWidth',2); title('x_{true}'); set(gca,'FontSize',20); 
hold off;
legend('short source','long source');
subplot(2,1,2); plot(s_vals,y_true,'-k','LineWidth',2); title('y_{true}'); set(gca,'FontSize',20);
hold on;
plot(s_vals,y_mod,'--m','LineWidth',2); title('y_{true}'); set(gca,'FontSize',20);
hold off;


% Compute ``naive'' solution using PREVIOUS model operator matrix G in 
% order to account for model/data misfit
x_solmod = G\y_datamod;

figure(5);
plot(t_mod,x_mod,'-k','LineWidth',2); hold on;
plot(t_mod,x_solmod,'--r','LineWidth',2); hold off
title('Naive solution using previous model G');
legend('x_{true}','x_{sol}');
set(gca,'FontSize',20);
xlim([0,1]);




