function [ F ] = fun( p,X )
% This function describes the camera projection
% from the world plane to the image plane.
% This function is used for the Levenberg-Marquardt method.
K = [p(1) p(2) p(3); 0 p(4) p(5); 0 0 1]; % Intrinsic parameter matrix K
w = [p(6) p(7) p(8)]';
t = [p(9) p(10) p(11)]';
Wx = [0,-w(3),w(2);w(3),0,-w(1);-w(2),w(1),0];
phi = sqrt(w'*w);
R = 1 + Wx*sin(phi)/phi + (Wx.^2)*(1-cos(phi))/phi;
P = K * [R t];
[M,N] = size(X);
F = zeros(M,N);
for k = 1:M
Y = P * [X(k,:),0,1]';
F(k,1) = Y(1)/Y(3);
F(k,2) = Y(2)/Y(3);
end
end
