clear all
close all
clc

n = 64;
[G, y, x] = phillips(n);
[U,S,V] = svd(G);
singular_values = diag(S);

tolerance = 10^-6;
iterations = 0;
error = 1;
X = randn(64,1);
errors = error;

while error > tolerance
    
    iterations = iterations + 1;
    X = G' * G * X / norm(X);
    rhos(iterations) = sqrt(norm(X));
    error = norm(rhos(iterations) - singular_values(1));
    errors(iterations) = error;
    
end
figure(1)
plot(1:length(rhos), rhos,'*', 1:length(rhos),...
singular_values(1)*(ones(length(rhos))))
figure(2)
plot(1:iterations, log(errors))