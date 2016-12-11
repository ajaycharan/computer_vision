clear; close;
%% Preprocess %%
% Read images
im1 = imread('WP_20160217_003.jpg');
im2 = imread('WP_20160217_004.jpg');
% Convert RGB image to grayscale
im1g = single(rgb2gray(im1));
im2g = single(rgb2gray(im2));
[R,C] = size(im1g);
%% SIFT of im1 and im2 by VL_SIFT function %%
[F1,D1] = vl_sift(im1g); % Each column of F is a feature frame
[F2,D2] = vl_sift(im2g); % Each column of D is a discriptor
d = dist(D1',D2); % Distance between D1's column and D2's column
[Y I] = min(d);
count = 0; % Number of non-overlapped correspondences
c1 = zeros(1,2); % Corresponding feature coordinates of im1
c2 = zeros(1,2); % Corresponding feature coordinates of im2
%% Find correspondences between two images using Euclidean distance %%
img = [im1,im2];
for k = 1:length(Y)
 ind = 1; % Indicator to avoid overlapped correspondences
 for l = 1:length(I)
 if l~=k && I(l)==I(k)
 ind = 0;
 break;
 end
 end
 if ind && Y(k) < 35 % Threshold for Euclidean distance
 count = count + 1;
 c1(count,:) = round(F1(1:2,I(k)));
 c2(count,:) = round(F2(1:2,k));
 end
end
%% RANSAC algorithm %%
nc = 6; % Number of correspondences used to find a homography
N = fix(log(1-.99)/log(1-(1-.1)^nc)); % Number of trials by 10% rule
M = fix((1-.1)*count); % Minimum size for the inlier set
d_min = 1e100;
for n = 1:N
 lcv = 1; % Loop control variable
 while lcv % To avoid repeated selection
 r = randi(count,nc,1);
 r = sort(r);
 for k = 1:nc-1
 lcv = lcv*(r(k+1)-r(k));
 end
 lcv = ~lcv;
 end
 A = zeros(2*nc,9);
 for k = 1:nc
 A(2*k-1:2*k,:)=...
 [0,0,0,-[c1(r(k),:),1],c2(r(k),2)*[c1(r(k),:),1];
 [c1(r(k),:),1],0,0,0,-c2(r(k),1)*[c1(r(k),:),1]];
 end
 [U,D,V] = svd(A);
 h = V(:,9);
 H = [h(1),h(2),h(3);h(4),h(5),h(6);h(7),h(8),h(9)];

 d2 = zeros(count,1); % d^2(x_measured, x_true)
 for k = 1:count
 x_true = H*[c1(k,:),1]'; % x_true in HC
 temp = x_true/x_true(3);
 x_true = temp(1:2); % x_true in image plane
 d = c2(k,:)-x_true';
 d2(k) = d(1)^2+d(2)^2;
 end
 [Y I] = sort(d2);
 if sum(Y(1:M)) < d_min
 d_min = sum(Y(1:M));
 inliers = I(1:M);
 outliers = I(M+1:end);
 end
end
% Visualize the inliers and outliers
figure; image(img); truesize; hold on;
for k = inliers'
 plot([c1(k,1),C+c2(k,1)],[c1(k,2),c2(k,2)],'-og','linewidth',1);
end
for k = outliers'
 plot([c1(k,1),C+c2(k,1)],[c1(k,2),c2(k,2)],'-or','linewidth',1);
end
plot([C,C],[1,R],'-k'); hold off;
%% Linear Least Squares %%
A = zeros(2*M,9);
for k = 1:M
 A(2*k-1:2*k,:)=...
 [0,0,0,-[c1(inliers(k),:),1],c2(inliers(k),2)*[c1(inliers(k),:),1];
 [c1(inliers(k),:),1],0,0,0,-c2(inliers(k),1)*[c1(inliers(k),:),1]];
end
[U,D,V] = svd(A);
h1 = V(:,9); % Homography estimated by LLS with all inliers
%% Non-linear Least Square (Levenberg-Marquardt) %%
c1 = c1(inliers,:)';
c1 = c1(:);
c2 = c2(inliers,:)';
c2 = c2(:);
opt = optimset('Algorithm','levenberg-marquardt');
h2 = lsqcurvefit(@fun,h1,c1,c2,[],[],opt); % Refined homography by L.M.
H = [h2(1),h2(2),h2(3);h2(4),h2(5),h2(6);h2(7),h2(8),h2(9)];
% Source Code for Image Mosaicking
% Mosaic function:

% Fun function:

Main Script:
clear; close; clc;
%% Preprocess %%
% Load single images.
im1 = imread('WP_20160217_003.jpg');
im2 = imread('WP_20160217_004.jpg');
im3 = imread('WP_20160217_005.jpg');
im4 = imread('WP_20160217_006.jpg');
im5 = imread('WP_20160217_007.jpg');
im6 = imread('WP_20160217_008.jpg');
im7 = imread('WP_20160217_009.jpg');
[M,N,C] = size(im2);
% Load estimated and refined homographies in previous steps.
% All the refined homographies were saved as mat files.
H12 = load('H12'); H12 = H12.H; % Homography of im1 to im2
H23 = load('H23'); H23 = H23.H; % Homography of im3 to im2
H34 = load('H34'); H34 = H34.H; % Homography of im1 to im2
H54 = load('H54'); H54 = H54.H; % Homography of im3 to im2
H65 = load('H65'); H65 = H65.H; % Homography of im1 to im2
H76 = load('H76'); H76 = H76.H; % Homography of im3 to im2
% im1 im2 im3 im4 im5 im6 im7 <- order of single images : im4 is in center.

H14 = H12*H23*H34;
H24 = H23*H34;
H64 = H65*H54;
H74 = H76*H65*H54;

%% Boundary Condition of Mosaiced Image %%
h14 = H14'; h14 = h14(:); % Change homograpy to a vector form.
h24 = H24'; h24 = h24(:);
h34 = H34'; h34 = h34(:);
h54 = H54'; h54 = h54(:);
h64 = H64'; h64 = h64(:);
h74 = H74'; h74 = h74(:);
c14 = fun(h14,[1,1,N,1,1,M,N,M]); % Transformed boundaries of im1
c74 = fun(h74,[1,1,N,1,1,M,N,M]); % Transformed boundaries of im7
x = [1,3,5,7];
y = [2,4,6,8];
xmin = round(min([c14(x);c74(x)]));
xmax = round(max([c14(x);c74(x)]));
ymin = round(min([c14(y);c74(y)]));
ymax = round(max([c14(y);c74(y)]));
%% Assign pixel values into the mosaiced image %%
img = zeros(ymax-ymin+1,xmax-xmin+1,C); % Initialize mosaiced image
img = mosaic(img,im1,H14,xmin,ymin); % Mosaicking im1
img = mosaic(img,im7,H74,xmin,ymin); % Mosaicking im7
img = mosaic(img,im2,H24,xmin,ymin); % Mosaicking im2
img = mosaic(img,im6,H64,xmin,ymin); % Mosaicking im6
img = mosaic(img,im3,H34,xmin,ymin); % Mosaicking im3
img = mosaic(img,im5,H54,xmin,ymin); % Mosaicking im5
img(2-ymin:M+1-ymin,2-xmin:N+1-xmin,:) = im4; % Mosaicking im4
figure; imshow(img); imwrite(img,'Mosaic_apt');