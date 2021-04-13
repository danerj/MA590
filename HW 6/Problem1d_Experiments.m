%% Increasing Blur

clear all
close all
clc

N = 64;
[G,y,x] = blur(N,10,10);

figure(1)
X = reshape(x,N,N);
imagesc(X)
title('The Original Image "X"')
colormap gray
axis image

figure(2)
Y = reshape(y,N,N);
imagesc(Y)
colormap gray
title('Result "Y" of blurring "X" through transformation "G"')
tic
X_CGLS = cgls(y,G,1000);
toc
figure(3)
imagesc(reshape(X_CGLS(:,11),N,N))
colormap gray
figure(4)
imagesc(reshape(X_CGLS(:,51),N,N))
colormap gray
figure(5)
imagesc(reshape(X_CGLS(:,101),N,N))
colormap gray
figure(6)
imagesc(reshape(X_CGLS(:,end),N,N))
colormap gray

%% Increasing Noise

clear all
close all
clc

N = 64;
[G,y,x] = blur(N,10,1.4);

figure(1)
X = reshape(x,N,N);
imagesc(X)
title('The Original Image "X"')
colormap gray
axis image

figure(2)
Y = reshape(y,N,N);
imagesc(Y)
colormap gray
title('Result "Y" of blurring "X" through transformation "G"')

noise = randn(length(y),1);
error = (0.5*norm(y))*noise/norm(noise); %Change 0.1 to other values to experiment with the noise level.

y_noisy = y + error;
n = 1000;
tic
X_CGLS_noisy = cgls(y_noisy, G, n);
toc
X_end_noisy = reshape(X_CGLS_noisy(:,end), N,N);

figure(3)
imagesc(X_end_noisy)
colormap gray

figure(4)
imagesc(reshape(X_CGLS_noisy(:,3), N,N))
colormap gray

figure(5)
imagesc(reshape(X_CGLS_noisy(:,3), N,N))
colormap gray

figure(6)
imagesc(reshape(X_CGLS_noisy(:,3), N,N))
colormap gray

residual_norms = zeros(1,n+2);
x_norms = zeros(1,n+2);

for ii = 1:n
    
    residual_norms(ii) = norm(G*X_CGLS_noisy(:,ii) - y);
    x_norms(ii) = norm(X_CGLS_noisy(:,ii));
    
end

loglog(residual_norms, x_norms);

x_errors = zeros(1,1002);
for kk = 1:1002
x_errors(kk) = norm(X_CGLS_noisy(:,kk) - x);
end
[minimum_error,minimum_error_index] = min(x_errors);
figure(7)
imagesc(reshape(X_CGLS_noisy(:,minimum_error_index),N,N))
colormap gray

%% Increasing Blur and Noise

clear all
close all
clc

N = 64;
[G,y,x] = blur(N,10,10);

figure(1)
X = reshape(x,N,N);
imagesc(X)
title('The Original Image "X"')
colormap gray
axis image

figure(2)
Y = reshape(y,N,N);
imagesc(Y)
colormap gray
title('Result "Y" of blurring "X" through transformation "G"')

noise = randn(length(y),1);
error = (1*norm(y))*noise/norm(noise); %Change 0.1 to other values to experiment with the noise level.

y_noisy = y + error;
n = 1000;
tic
X_CGLS_noisy = cgls(y_noisy, G, n);
toc
X_end_noisy = reshape(X_CGLS_noisy(:,end), N,N);

figure(3)
imagesc(X_end_noisy)
colormap gray

figure(4)
imagesc(reshape(X_CGLS_noisy(:,3), N,N))
colormap gray

figure(5)
imagesc(reshape(X_CGLS_noisy(:,3), N,N))
colormap gray

figure(6)
imagesc(reshape(X_CGLS_noisy(:,3), N,N))
colormap gray

residual_norms = zeros(1,n+2);
x_norms = zeros(1,n+2);

for ii = 1:n
    
    residual_norms(ii) = norm(G*X_CGLS_noisy(:,ii) - y);
    x_norms(ii) = norm(X_CGLS_noisy(:,ii));
    
end

loglog(residual_norms, x_norms);

x_errors = zeros(1,1002);
for kk = 1:1002
x_errors(kk) = norm(X_CGLS_noisy(:,kk) - x);
end
[minimum_error,minimum_error_index] = min(x_errors);
figure(7)
imagesc(reshape(X_CGLS_noisy(:,minimum_error_index),N,N))
colormap gray
