clc
clear all
close all
% toatal no. of mugshots
N = 60;
% width of image
w = 92;
% height of image
h = 112;
% size of image
sz = w*h;
% matrix A containig column vectors as image vectors
imvec = zeros(sz, N);
% read the images and convert them to columns and store in A
% imshow(imread(('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s1\1.pgm')));
imvec(:,1) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s1\1.pgm');
imvec(:,2) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s1\2.pgm');
imvec(:,3) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s1\3.pgm');
imvec(:,4) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s1\4.pgm');
imvec(:,5) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s1\5.pgm');
imvec(:,6) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s1\6.pgm');
imvec(:,7) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s1\7.pgm');
imvec(:,8) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s1\8.pgm');
imvec(:,9) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s1\9.pgm');
imvec(:,10) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s1\10.pgm');
imvec(:,11) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s2\1.pgm');
imvec(:,12) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s2\2.pgm');
imvec(:,13) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s2\3.pgm');
imvec(:,14) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s2\4.pgm');
imvec(:,15) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s2\5.pgm');
imvec(:,16) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s2\6.pgm');
imvec(:,17) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s2\7.pgm');
imvec(:,18) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s2\8.pgm');
imvec(:,19) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s2\9.pgm');
imvec(:,20) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s2\10.pgm');
imvec(:,21) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s3\1.pgm');
imvec(:,22) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s3\2.pgm');
imvec(:,23) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s3\3.pgm');
imvec(:,24) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s3\4.pgm');
imvec(:,25) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s3\5.pgm');
imvec(:,26) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s3\6.pgm');
imvec(:,27) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s3\7.pgm');
imvec(:,28) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s3\8.pgm');
imvec(:,29) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s3\9.pgm');
imvec(:,30) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s3\10.pgm');
imvec(:,31) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s4\1.pgm');
imvec(:,32) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s4\2.pgm');
imvec(:,33) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s4\3.pgm');
imvec(:,34) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s4\4.pgm');
imvec(:,35) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s4\5.pgm');
imvec(:,36) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s4\6.pgm');
imvec(:,37) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s4\7.pgm');
imvec(:,38) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s4\8.pgm');
imvec(:,39) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s4\9.pgm');
imvec(:,40) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s4\10.pgm');
imvec(:,41) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s5\1.pgm');
imvec(:,42) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s5\2.pgm');
imvec(:,43) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s5\3.pgm');
imvec(:,44) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s5\4.pgm');
imvec(:,45) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s5\5.pgm');
imvec(:,46) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s5\6.pgm');
imvec(:,47) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s5\7.pgm');
imvec(:,48) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s5\8.pgm');
imvec(:,49) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s5\9.pgm');
imvec(:,50) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s5\10.pgm');
imvec(:,51) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s6\1.pgm');
imvec(:,52) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s6\2.pgm');
imvec(:,53) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s6\3.pgm');
imvec(:,54) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s6\4.pgm');
imvec(:,55) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s6\5.pgm');
imvec(:,56) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s6\6.pgm');
imvec(:,57) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s6\7.pgm');
imvec(:,58) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s6\8.pgm');
imvec(:,59) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s6\9.pgm');
imvec(:,60) = imVecNg('C:\Users\Ajay Charan\Desktop\ASAP\att_faces\orl_faces\s6\10.pgm');

% matrix to get correlation
C = imvec'*imvec;

% eigen decomposition
[vectorC, valueC] = eig(C);

% eigen value matrix with diagonal elements as eigen values
ss = sum(valueC);

% sort the eigen value matrix
[ss, iii] = sort(ss, 'descend');
sizI = size(iii);

% sort the corresponding eigen vectors
for i = 1:10
    
    vectorCInt(:,i) = vectorC(:, iii(i));
end

vectorC = vectorCInt;

% choose the most significant eigen vectors
vectorL = imvec*vectorC;
% vectorL = vectorL(1:10);

% 
Coeff = imvec'*vectorL;

for i = 1:6
    model(i,:) = mean(Coeff((10*(i-1)+1):10*i,:));
end

while(1)
    imagename = input('Enter the filename of image to recognize: ','s');
    switch imagename(1)
        case '1'
            imagenameNum = str2num(imagename(2));
        case '2'
            imagenameNum = 10+str2num(imagename(2));
        case '3'
            imagenameNum = 20+str2num(imagename(2));
        case '4'
            imagenameNum = 30+str2num(imagename(2));
        case '5'
            imagenameNum = 40+str2num(imagename(2));
        case '6'
            imagenameNum = 50+str2num(imagename(2));
    end
    imageco = imvec(:,imagenameNum)'*vectorL;
    imshow(uint8(vecIm(imvec(:,imagenameNum),w,h)));
    top = 1;
    for i = 1:6
        if(norm(model(i,:)-imageco,1)< norm(model(top,:)-imageco,1))
            top = i;
        end
    end
    message = sprintf('the image was of person %d', top);
    disp(message);
end
