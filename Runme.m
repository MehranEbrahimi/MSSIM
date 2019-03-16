%--------------------------------------------------------------------------
% Note: don't click on other figures while doing the iterations because
% iterations are displayed in the current figure.
%--------------------------------------------------------------------------
close all;
A = (imread('cameraman.tif'));


A_noisy = (imnoise(A, 'gaussian', 0, 1*0.005));


[mssim, mssim_derivative]=ssimwithderivative(double(A),double(A_noisy));
Q=double(A);
fctn = @(x) ssimwithderivative(reshape(double(x),[256,256]),double(A_noisy));
[fig,ph,th] = checkDerivative(fctn,Q(:));


lambda_1 = 0.01;
lambda_2 = 100; %
figure;
subplot(411); imshow(A); title('Original');
subplot(412); imshow(A_noisy); title('Noisy input');
subplot(413); u1=totalvar(A_noisy, lambda_1); title(['lambda = ', num2str(lambda_1)]); 
subplot(414); u2=totalvar(A_noisy, lambda_2); title(['lambda = ', num2str(lambda_2)]);

ssim(A_noisy,A)
ssim(u1,A)
ssim(u2,A)

psnr(A_noisy,A)
psnr(u1,A)
psnr(u2,A)
