clear all
close all
clc

%Setting up the Wing model test problem.
n = 100;
[G, y_true, x_true] = wing(n);

%Collecting SVD of G.
[U,S,V] = svd(G);
singular_values = diag(S);


for k = 1:n

   %Calculating the TSVD solution for each truncation parameter k.
   x_k =  V(:,1:k)*inv(S(1:k,1:k))*(U(:,1:k))'*y_true;
   plot(1:n, x_k, 1:n, x_true)
   title({'k = ' k})
   pause(1)
   
end

