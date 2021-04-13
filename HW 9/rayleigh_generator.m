function [ rayleigh_values ] = rayleigh_generator( sigma, n)
%rayleigh_generator.m produces a specified number of points following a
%Rayleight Distribution. This function omits the intermediate step of
%creating a vector valued random variable A = X*e_1 + Y*e_2 since this step
%is unnecessary in producing values from the distribution and causes a
%large loss in efficiency. (Unnecessary because ||A||_2 can be calculated
%using the first two entries of A, which are X and Y and are the only
%elements of A that are not guaranteed to be zero by definition. Therefore,
%Rayleigh values, R, are of the form R = sqrt(X^2 + Y^2) for X,Y normally
%distributed random variables with mean 0 and variance sigma^2.

%   Detailed explanation goes here

    %sigma - standard deviation of normal distribution from which we draw
    %to distinct points, X and Y, for producing each individual point from
    %the Rayleigh Distribution. 
    
    %n - the number of points from the Rayleigh Distribution to produce. 

rayleigh_values = zeros(n,1);

for ii = 1:n
    
    X = sigma*randn(1);
    Y = sigma*randn(1);
    rayleigh_values(ii) = sqrt(X^2+Y^2);

end

