clear all
close all
clc

N = 64;
[G,y,x] = blur(N);

noise = randn(length(y),1);
error = (0.1*norm(y))*noise/norm(noise); %Change 0.1 to other values to experiment with the noise level.

y_noise  = y + error;

GTG = G'*G;
GTy = G'*y_noise;

X = zeros(length(x),1);

change = 1;
tolerance = 10^-2;
iteration = 1;

while change > tolerance
    
    p_k = GTG*X(:,iteration) - GTy;
    alpha_k = (norm(p_k))^2/(norm(G*p_k))^2;
    X(:,iteration+1) = X(:,iteration) - alpha_k*p_k;
    change = norm(X(:,iteration+1) - X(:,iteration));
    iteration = iteration + 1;
    
end

