clear all; close all; clc;
imNum = '010203040506070810111213141516171819204748495051525354555657585960166263646566';
HOMO = cell(1,length(imNum)/2); % Cell to store homographies
CON = cell(1,length(imNum)/2); % Cell to store the coordinates of corners
FNAME = cell(1,length(imNum)/2); % Cell to store file names
for Num = 1:length(imNum)/2;
%% Image Load %%
fname = ['P10100',imNum(2*Num-1:2*Num),'s.jpg'];
FNAME{Num} = fname;
im_RGB = imread(fname); % RGB color image
im_gray = rgb2gray(im_RGB); % Gray-scale image
[M,N] = size(im_gray); % Size of image
%% Canny Edge Detection and Hough Transform %%
BW = edge(im_gray,'canny',0.7); % Binary image from Canny detector
% figure; imshow(BW);
[H,T,R] = hough(BW); % Hough transform
P = houghpeaks(H,18,'Threshold',1);
lines = houghlines(BW,T,R,P,'FillGap',100,'MinLength',100);
% Plot Hough lines
x = 1:N;
ab = zeros(length(lines),2);
% figure; imshow(im_RGB); hold on;
for k = 1:length(lines)
delta = lines(k).point2 - lines(k).point1;
a = delta(2)/delta(1); % Slope of line
if a == inf
% plot(lines(k).point1(1)*ones(1,M),1:M,'g');
b = inf;
else
b = lines(k).point1(2) - a*lines(k).point1(1);
% plot(x,a*x+b,'g');
end
ab(k,1) = a; % Save line parameters a and b
ab(k,2) = b; % where y = a*x + b
end
% hold off;
%% Corner Detection and Label %%
[S,I] = sort(abs(ab(:,1))); % Find horizontal lines
a_H = ab(I(1:10),1); % a of horizontal lines
b_H = ab(I(1:10),2); % b of horizontal lines
[b_H,I] = sort(b_H);
a_H = a_H(I);
Ico = corner(im_gray,80); % Corners in image plane
label = zeros(80,1); % Initialize integer label
for k = 1:10
d2 = (a_H(k)*Ico(:,1)-Ico(:,2)+b_H(k)).^2/(a_H(k)^2+1);
[d2,I] = sort(d2);
I = I(1:8);
[S,J] = sort(Ico(I,1));
label(8*(k-1)+1:8*k) = I(J);
end
CON{Num} = Ico(label,:);
% figure; imshow(im_RGB); hold on;
% for k = 1:80
% plot(Ico(label(k),1),Ico(label(k),2),'*r');
% text(Ico(label(k),1)+3,Ico(label(k),2),num2str(k),'color','y');
% end
% hold off;
%% Find Homograhy by LLS %%
for k = 1:10
Wco(8*(k-1)+1:8*k,1) = 0:24:168; % Corners in world plane (mm)
Wco(8*(k-1)+1:8*k,2) = 24*(k-1);
end
A = zeros(2*80,9);
for k = 1:80
A(2*k-1:2*k,:)=...
[0,0,0,-[Wco(k,:),1],Ico(label(k),2)*[Wco(k,:),1];
[Wco(k,:),1],0,0,0,-Ico(label(k),1)*[Wco(k,:),1]];
end
[U,D,V] = svd(A);
h = V(:,9); % Homography estimated by LLS with 80 correspondences
HOMO{Num} = vec2mat(h,3); % Save homographies
v12 = [h(1)*h(2),h(1)*h(5)+h(4)*h(2),h(4)*h(5),h(7)*h(2)+h(1)*h(8),...
h(7)*h(5)+h(4)*h(8),h(7)*h(8)];
v11 = [h(1)*h(1),h(1)*h(4)+h(4)*h(1),h(4)*h(4),h(7)*h(1)+h(1)*h(7),...
h(7)*h(4)+h(4)*h(7),h(7)*h(7)];
v22 = [h(2)*h(2),h(2)*h(5)+h(5)*h(2),h(5)*h(5),h(8)*h(2)+h(2)*h(8),...
h(8)*h(5)+h(5)*h(8),h(8)*h(8)];
B(2*Num-1:2*Num,:) = [v12;v11-v22];
end
[U,D,V] = svd(B);
v = V(:,6);
B = [v(1),v(2),v(4);v(2),v(3),v(5);v(4),v(5),v(6)];
%% Compute Intrinsic Parameters in K %%
y0 = (v(2)*v(4)-v(1)*v(5))/(v(1)*v(3)-v(2)^2); % Principle point y
lambda = v(6)-(v(4)^2+y0*(v(2)*v(4)-v(1)*v(5)))/v(1);
a_x = sqrt(lambda/v(1)); % Scale factor in x
a_y = sqrt(lambda*v(1)/(v(1)*v(3)-v(2)^2)); % Scale factor in y
s = -v(2)*a_x^2*a_y/lambda; % Skewness factor
x0 = s*y0/a_y-v(4)*a_x^2/lambda; % Principle point x
K = [a_x s x0; 0 a_y y0; 0 0 1]; % Intrinsic parameter matrix K
%% Compute Extrinsic Parameters %%
Rt = cell(1,Num); % Cell to store extrinsic parameters for each image
P = cell(1,Num); % Cell to store camera matrix for each image
for k = 1:Num
H = HOMO{k};
r1 = lambda*K\H(:,1);
r2 = lambda*K\H(:,2);
r3 = cross(r1,r2);
t = lambda*K\H(:,3);
R = [r1,r2,r3];
% phi = acos((trace(R)-1)/2);
% w = phi*[R(3,2)-R(2,3);R(1,3)-R(3,1);R(2,1)-R(1,2)]/(2*sin(phi));
% % Optimization by Levenberg-Marquardt method
% opt = optimset('Algorithm','levenberg-marquardt');
% p = [a_x,s,x0,a_y,y0,w(1),w(2),w(3),t(1),t(2),t(3)];
% q = lsqcurvefit(@fun,p,Wco,CON{k},[],[],opt); % Refined parameters
% K = [q(1) q(2) q(3); 0 q(4) q(5); 0 0 1];
% R = exp([0,-q(8),q(7);q(8),0,-q(6);-q(7),q(6),0]);
% t = [q(9),q(10),q(11)]';
Rt{k} = [R,t]; % Extrinsic parameter matrix [R|t]
P{k} = K*Rt{k}
end
%% Test Accuracy by Re-projecting Corner Points %%
n1 = 1; % Image Number having corners to be reprojected
n2 = 3; % Image Number which the corners reprojected to
im = imread(FNAME{n2});
figure; imshow(im); hold on;
for k = 1:80
plot(CON{n2}(k,1),CON{n2}(k,2),'*r');
X = P{n1}'/(P{n1}*P{n1}')*[CON{n1}(k,:),1]';
x = P{n2}*X;
plot(x(1)/x(3),x(2)/x(3),'*g');
end
hold off;
