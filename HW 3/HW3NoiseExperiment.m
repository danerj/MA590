%%
%Noise std dev 10^-3

for ii=1:1
    
clear all
close all
clc
%a

n = 64;

%The phillips.m function produces Phillips’ test problem, which is a
%discretization of a Fredholm integral equation of the first kind.
[G, y_true, x_true] = phillips(64);

%Note that G is poorly conditioned.
G_condition = cond(G);

%We generate noisy data y.
e = 10^-3*randn(length(y_true),1);
norm_e_test = norm(e);
y = y_true + e; %Noisy data

%Note that the noise does not appear at first to make a large difference.
plot1 = plot(x_true, y_true, '+', x_true, y, 'o');
xlabel('$$x_{true}$$','FontSize',16,'interpreter','latex');
ylabel('$$y, y_{true}$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 1: Initial Comparison between True and Noisy Solution' ''},'Interpreter','latex')
names = {'x_{true} vs. y_{true}' , 'x_{true} vs. y'};
legend(plot1,names, 'location', 'northwest')

%b

%The full svd of G.
[U,S,V] = svd(G);

%Analyze the singular values of G.
%Note that all singular values are nonzero, although many are 'small'.
singular_values = diag(S);
figure(2)
plot2 = semilogy(1:1:length(singular_values), singular_values, 'g*');
xlabel('$$i$$','FontSize',16,'interpreter','latex');
ylabel('$$\sigma_i$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 2: A Visual Inspection of the Singular Values of G' ''},'Interpreter','latex')
largest_singular_value = max(singular_values);
smallest_singular_value = min(singular_values);


%c

coefficients = abs(U'*y);
ratios = coefficients ./ singular_values;

figure(3)
plot3 = semilogy(1:64, coefficients,':', 1:64, singular_values, 'g*', 1:64, ratios, 'c-.');
xlabel('index $$i$$','FontSize',16,'interpreter','latex');
ylabel('coefficients, singular values, ratios','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 3: Graphically determining whether the Discrete Picard Condition holds' ''},'Interpreter','latex')
names = {'|u_i^Ty|' , '\sigma_i', '|u_i^Ty|/\sigma_i'};
legend(plot3,names, 'location', 'north')
%DPC Holds Here? Appears that yes it does since the coefficients tend to
%grow faster than the singular values as the index increases up to n=64.

%d

%There are no nonzero eigenvalues in our case. So p = n. Still we write the
%code to be flexible in case changing earlier steps results in zero
%singular values in future experiments. This means that the compact svd is
%the same as the full svd.
p = n;
nonzero_singular_values = singular_values(1:p);

U_c = U(:,1:p);
S_c = S(1:p,1:p);
V_c = V(:,1:p);

%The pseudoinverse of G.
G_dagger = V_c*inv(S_c)*(U_c)';

%The generalized inverse solution.
x_dagger = G_dagger*y;

figure(4)
plot4 = plot(1:n, x_true, 'm--', 1:n, x_dagger, 'k:');
xlabel('index $$i$$','FontSize',16,'interpreter','latex');
ylabel('$$x_{true}, x^+$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 4: Considering the limited practicality of the generalized inverse solution' ''},'Interpreter','latex')
names = {'x_{true}' , 'x^+'};
legend(plot4,names, 'location', 'southwest');


%e

x_k_norms = zeros(1,n);
residual_norms = zeros(1,n);

%Finding approximate solutions using truncated svd. It is simpler to to
%calculate x_k for k=1,...,64 even though we will only inspect k=2,..62.
for ii = 1:n

   x_k =  V(:,1:ii)*inv(S(1:ii,1:ii))*(U(:,1:ii))'*y;
   x_k_norms(ii) = norm(x_k);
   residual_norms(ii) = norm(G*x_k - y);
   
end

%We see that the 'corner' of the L-curve is at about the point where our
%vertical axis value is about 3. The 2-norm of x_k first surpasses 3 at
%about k = 7. The L-Curve method is not an exact calculation but a
%heuristic method, so the truly optimal parameter may be slightly more or
%less than 7.
figure(5)
plot5 = loglog(residual_norms(2:62), x_k_norms(2:62), 'm:');
xlabel('$$||e_k||_2$$','FontSize',16,'interpreter','latex');
ylabel('$$||x_k||_2$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 5: Finding optimal truncation parameter with an L-curve' ''},'Interpreter','latex')

% solution using truncation param k = 7

k_optimal = 7;

x_k_optimal =  V(:,1:k_optimal)*inv(S(1:k_optimal,1:k_optimal))*(U(:,1:k_optimal))'*y;
figure(6)
plot6 = plot(1:n, x_true, 'm--', 1:n, x_k_optimal, 'g:');
xlabel('index $$i$$','FontSize',16,'interpreter','latex');
ylabel('$$x_{true}, x_{k_{optimal}}$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 6: Comparing true solution and optimally truncated approximate solution' ''},'Interpreter','latex')
names = {'x_{true}' , 'x_{k_{optimal}}'};
legend(plot6,names, 'location', 'northeast')

end

%%
%Noise std dev 10^-2
for ii=1:1

clear all
close all
clc

%a

n = 64;

%The phillips.m function produces Phillips’ test problem, which is a
%discretization of a Fredholm integral equation of the first kind.
[G, y_true, x_true] = phillips(64);

%Note that G is poorly conditioned.
G_condition = cond(G);

%We generate noisy data y.
e = 10^-2*randn(length(y_true),1);
norm_e_test = norm(e);
y = y_true + e; %Noisy data

%Note that the noise does not appear at first to make a large difference.
plot1 = plot(x_true, y_true, '+', x_true, y, 'o');
xlabel('$$x_{true}$$','FontSize',16,'interpreter','latex');
ylabel('$$y, y_{true}$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 1: Initial Comparison between True and Noisy Solution' ''},'Interpreter','latex')
names = {'x_{true} vs. y_{true}' , 'x_{true} vs. y'};
legend(plot1,names, 'location', 'northwest')

%b

%The full svd of G.
[U,S,V] = svd(G);

%Analyze the singular values of G.
%Note that all singular values are nonzero, although many are 'small'.
singular_values = diag(S);
figure(2)
plot2 = semilogy(1:1:length(singular_values), singular_values, 'g*');
xlabel('$$i$$','FontSize',16,'interpreter','latex');
ylabel('$$\sigma_i$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 2: A Visual Inspection of the Singular Values of G' ''},'Interpreter','latex')
largest_singular_value = max(singular_values);
smallest_singular_value = min(singular_values);


%c

coefficients = abs(U'*y);
ratios = coefficients ./ singular_values;

figure(3)
plot3 = semilogy(1:64, coefficients,':', 1:64, singular_values, 'g*', 1:64, ratios, 'c-.');
xlabel('index $$i$$','FontSize',16,'interpreter','latex');
ylabel('coefficients, singular values, ratios','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 3: Graphically determining whether the Discrete Picard Condition holds' ''},'Interpreter','latex')
names = {'|u_i^Ty|' , '\sigma_i', '|u_i^Ty|/\sigma_i'};
legend(plot3,names, 'location', 'north')
%DPC Holds Here? Appears that yes it does since the coefficients tend to
%grow faster than the singular values as the index increases up to n=64.

%d

%There are no nonzero eigenvalues in our case. So p = n. Still we write the
%code to be flexible in case changing earlier steps results in zero
%singular values in future experiments. This means that the compact svd is
%the same as the full svd.
p = n;
nonzero_singular_values = singular_values(1:p);

U_c = U(:,1:p);
S_c = S(1:p,1:p);
V_c = V(:,1:p);

%The pseudoinverse of G.
G_dagger = V_c*inv(S_c)*(U_c)';

%The generalized inverse solution.
x_dagger = G_dagger*y;

figure(4)
plot4 = plot(1:n, x_true, 'm--', 1:n, x_dagger, 'k:');
xlabel('index $$i$$','FontSize',16,'interpreter','latex');
ylabel('$$x_{true}, x^+$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 4: Considering the limited practicality of the generalized inverse solution' ''},'Interpreter','latex')
names = {'x_{true}' , 'x^+'};
legend(plot4,names, 'location', 'southwest');


%e

x_k_norms = zeros(1,n);
residual_norms = zeros(1,n);

%Finding approximate solutions using truncated svd. It is simpler to to
%calculate x_k for k=1,...,64 even though we will only inspect k=2,..62.
for ii = 1:n

   x_k =  V(:,1:ii)*inv(S(1:ii,1:ii))*(U(:,1:ii))'*y;
   x_k_norms(ii) = norm(x_k);
   residual_norms(ii) = norm(G*x_k - y);
   
end

%We see that the 'corner' of the L-curve is at about the point where our
%vertical axis value is about 3. The 2-norm of x_k first surpasses 3 at
%about k = 7. The L-Curve method is not an exact calculation but a
%heuristic method, so the truly optimal parameter may be slightly more or
%less than 7.
figure(5)
plot5 = loglog(residual_norms(2:62), x_k_norms(2:62), 'm:');
xlabel('$$||e_k||_2$$','FontSize',16,'interpreter','latex');
ylabel('$$||x_k||_2$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 5: Finding optimal truncation parameter with an L-curve' ''},'Interpreter','latex')

% solution using truncation param k = 7

k_optimal = 7;

x_k_optimal =  V(:,1:k_optimal)*inv(S(1:k_optimal,1:k_optimal))*(U(:,1:k_optimal))'*y;
figure(6)
plot6 = plot(1:n, x_true, 'm--', 1:n, x_k_optimal, 'g:');
xlabel('index $$i$$','FontSize',16,'interpreter','latex');
ylabel('$$x_{true}, x_{k_{optimal}}$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 6: Comparing true solution and optimally truncated approximate solution' ''},'Interpreter','latex')
names = {'x_{true}' , 'x_{k_{optimal}}'};
legend(plot6,names, 'location', 'northeast')

end

%%
%Noise std dev 10^-1

for ii=1:1
    
clear all
close all
clc

%a

n = 64;

%The phillips.m function produces Phillips’ test problem, which is a
%discretization of a Fredholm integral equation of the first kind.
[G, y_true, x_true] = phillips(64);

%Note that G is poorly conditioned.
G_condition = cond(G);

%We generate noisy data y.
e = 10^-1*randn(length(y_true),1);
norm_e_test = norm(e);
y = y_true + e; %Noisy data

%Note that the noise does not appear at first to make a large difference.
plot1 = plot(x_true, y_true, '+', x_true, y, 'o');
xlabel('$$x_{true}$$','FontSize',16,'interpreter','latex');
ylabel('$$y, y_{true}$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 1: Initial Comparison between True and Noisy Solution' ''},'Interpreter','latex')
names = {'x_{true} vs. y_{true}' , 'x_{true} vs. y'};
legend(plot1,names, 'location', 'northwest')

%b

%The full svd of G.
[U,S,V] = svd(G);

%Analyze the singular values of G.
%Note that all singular values are nonzero, although many are 'small'.
singular_values = diag(S);
figure(2)
plot2 = semilogy(1:1:length(singular_values), singular_values, 'g*');
xlabel('$$i$$','FontSize',16,'interpreter','latex');
ylabel('$$\sigma_i$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 2: A Visual Inspection of the Singular Values of G' ''},'Interpreter','latex')
largest_singular_value = max(singular_values);
smallest_singular_value = min(singular_values);


%c

coefficients = abs(U'*y);
ratios = coefficients ./ singular_values;

figure(3)
plot3 = semilogy(1:64, coefficients,':', 1:64, singular_values, 'g*', 1:64, ratios, 'c-.');
xlabel('index $$i$$','FontSize',16,'interpreter','latex');
ylabel('coefficients, singular values, ratios','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 3: Graphically determining whether the Discrete Picard Condition holds' ''},'Interpreter','latex')
names = {'|u_i^Ty|' , '\sigma_i', '|u_i^Ty|/\sigma_i'};
legend(plot3,names, 'location', 'north')
%DPC Holds Here? Appears that yes it does since the coefficients tend to
%grow faster than the singular values as the index increases up to n=64.

%d

%There are no nonzero eigenvalues in our case. So p = n. Still we write the
%code to be flexible in case changing earlier steps results in zero
%singular values in future experiments. This means that the compact svd is
%the same as the full svd.
p = n;
nonzero_singular_values = singular_values(1:p);

U_c = U(:,1:p);
S_c = S(1:p,1:p);
V_c = V(:,1:p);

%The pseudoinverse of G.
G_dagger = V_c*inv(S_c)*(U_c)';

%The generalized inverse solution.
x_dagger = G_dagger*y;

figure(4)
plot4 = plot(1:n, x_true, 'm--', 1:n, x_dagger, 'k:');
xlabel('index $$i$$','FontSize',16,'interpreter','latex');
ylabel('$$x_{true}, x^+$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 4: Considering the limited practicality of the generalized inverse solution' ''},'Interpreter','latex')
names = {'x_{true}' , 'x^+'};
legend(plot4,names, 'location', 'southwest');


%e

x_k_norms = zeros(1,n);
residual_norms = zeros(1,n);

%Finding approximate solutions using truncated svd. It is simpler to to
%calculate x_k for k=1,...,64 even though we will only inspect k=2,..62.
for ii = 1:n

   x_k =  V(:,1:ii)*inv(S(1:ii,1:ii))*(U(:,1:ii))'*y;
   x_k_norms(ii) = norm(x_k);
   residual_norms(ii) = norm(G*x_k - y);
   
end

%We see that the 'corner' of the L-curve is at about the point where our
%vertical axis value is about 3. The 2-norm of x_k first surpasses 3 at
%about k = 7. The L-Curve method is not an exact calculation but a
%heuristic method, so the truly optimal parameter may be slightly more or
%less than 7.
figure(5)
plot5 = loglog(residual_norms(2:62), x_k_norms(2:62), 'm:');
xlabel('$$||e_k||_2$$','FontSize',16,'interpreter','latex');
ylabel('$$||x_k||_2$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 5: Finding optimal truncation parameter with an L-curve' ''},'Interpreter','latex')

% solution using truncation param k = 7

k_optimal = 7;

x_k_optimal =  V(:,1:k_optimal)*inv(S(1:k_optimal,1:k_optimal))*(U(:,1:k_optimal))'*y;
figure(6)
plot6 = plot(1:n, x_true, 'm--', 1:n, x_k_optimal, 'g:');
xlabel('index $$i$$','FontSize',16,'interpreter','latex');
ylabel('$$x_{true}, x_{k_{optimal}}$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 6: Comparing true solution and optimally truncated approximate solution' ''},'Interpreter','latex')
names = {'x_{true}' , 'x_{k_{optimal}}'};
legend(plot6,names, 'location', 'northeast')

end

%%

for ii=1:1
    
clear all
close all
clc
%a

n = 64;

%The phillips.m function produces Phillips’ test problem, which is a
%discretization of a Fredholm integral equation of the first kind.
[G, y_true, x_true] = phillips(64);

%Note that G is poorly conditioned.
G_condition = cond(G);

%We generate noisy data y.
e = 10^0*randn(length(y_true),1);
norm_e_test = norm(e);
y = y_true + e; %Noisy data

%Note that the noise does not appear at first to make a large difference.
plot1 = plot(x_true, y_true, '+', x_true, y, 'o');
xlabel('$$x_{true}$$','FontSize',16,'interpreter','latex');
ylabel('$$y, y_{true}$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 1: Initial Comparison between True and Noisy Solution' ''},'Interpreter','latex')
names = {'x_{true} vs. y_{true}' , 'x_{true} vs. y'};
legend(plot1,names, 'location', 'northwest')

%b

%The full svd of G.
[U,S,V] = svd(G);

%Analyze the singular values of G.
%Note that all singular values are nonzero, although many are 'small'.
singular_values = diag(S);
figure(2)
plot2 = semilogy(1:1:length(singular_values), singular_values, 'g*');
xlabel('$$i$$','FontSize',16,'interpreter','latex');
ylabel('$$\sigma_i$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 2: A Visual Inspection of the Singular Values of G' ''},'Interpreter','latex')
largest_singular_value = max(singular_values);
smallest_singular_value = min(singular_values);


%c

coefficients = abs(U'*y);
ratios = coefficients ./ singular_values;

figure(3)
plot3 = semilogy(1:64, coefficients,':', 1:64, singular_values, 'g*', 1:64, ratios, 'c-.');
xlabel('index $$i$$','FontSize',16,'interpreter','latex');
ylabel('coefficients, singular values, ratios','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 3: Graphically determining whether the Discrete Picard Condition holds' ''},'Interpreter','latex')
names = {'|u_i^Ty|' , '\sigma_i', '|u_i^Ty|/\sigma_i'};
legend(plot3,names, 'location', 'north')
%DPC Holds Here? Appears that yes it does since the coefficients tend to
%grow faster than the singular values as the index increases up to n=64.

%d

%There are no nonzero eigenvalues in our case. So p = n. Still we write the
%code to be flexible in case changing earlier steps results in zero
%singular values in future experiments. This means that the compact svd is
%the same as the full svd.
p = n;
nonzero_singular_values = singular_values(1:p);

U_c = U(:,1:p);
S_c = S(1:p,1:p);
V_c = V(:,1:p);

%The pseudoinverse of G.
G_dagger = V_c*inv(S_c)*(U_c)';

%The generalized inverse solution.
x_dagger = G_dagger*y;

figure(4)
plot4 = plot(1:n, x_true, 'm--', 1:n, x_dagger, 'k:');
xlabel('index $$i$$','FontSize',16,'interpreter','latex');
ylabel('$$x_{true}, x^+$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 4: Considering the limited practicality of the generalized inverse solution' ''},'Interpreter','latex')
names = {'x_{true}' , 'x^+'};
legend(plot4,names, 'location', 'southwest');


%e

x_k_norms = zeros(1,n);
residual_norms = zeros(1,n);

%Finding approximate solutions using truncated svd. It is simpler to to
%calculate x_k for k=1,...,64 even though we will only inspect k=2,..62.
for ii = 1:n

   x_k =  V(:,1:ii)*inv(S(1:ii,1:ii))*(U(:,1:ii))'*y;
   x_k_norms(ii) = norm(x_k);
   residual_norms(ii) = norm(G*x_k - y);
   
end

%We see that the 'corner' of the L-curve is at about the point where our
%vertical axis value is about 3. The 2-norm of x_k first surpasses 3 at
%about k = 7. The L-Curve method is not an exact calculation but a
%heuristic method, so the truly optimal parameter may be slightly more or
%less than 7.
figure(5)
plot5 = loglog(residual_norms(2:62), x_k_norms(2:62), 'm:');
xlabel('$$||e_k||_2$$','FontSize',16,'interpreter','latex');
ylabel('$$||x_k||_2$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 5: Finding optimal truncation parameter with an L-curve' ''},'Interpreter','latex')

% solution using truncation param k = 7

k_optimal = 7;

x_k_optimal =  V(:,1:k_optimal)*inv(S(1:k_optimal,1:k_optimal))*(U(:,1:k_optimal))'*y;
figure(6)
plot6 = plot(1:n, x_true, 'm--', 1:n, x_k_optimal, 'g:');
xlabel('index $$i$$','FontSize',16,'interpreter','latex');
ylabel('$$x_{true}, x_{k_{optimal}}$$','FontSize',16,'interpreter','latex'); % make the label look pretty
title({'' 'Figure 6: Comparing true solution and optimally truncated approximate solution' ''},'Interpreter','latex')
names = {'x_{true}' , 'x_{k_{optimal}}'};
legend(plot6,names, 'location', 'northeast')

end