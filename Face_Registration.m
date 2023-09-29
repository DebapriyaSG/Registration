%This code registers two images using the Mutual_Information_Computation
%function

clear;
clc; 

%Initial guess
x0=[-16,-251,22.5];
x0(1)=(pi/180)*x0(1); 

%Global optimizer

opt=optimoptions('fmincon','Display','iter','StepTolerance',1e-2);
gs=GlobalSearch;
problem=createOptimProblem('fmincon','x0',x0,'objective',@Mutual_Information_Computation,'lb',[(pi/180)*(-16.5),-253,22],'ub',[(pi/180)*(-14),-249,25],'options',opt);
[x,fval]=run(gs,problem);

%Final output

%Hard coded
A=imread('Data\Face1.jpg'); %Reference Image
B=imread('Data\Face1_Rotate_15_X250_Y-24.jpg'); %Moving Image

%Fetching data from location

% path1=strcat(pwd,'\Data');
% listing=dir(path1);
% reference=strcat(path1,'\',listing(3).name);
% moving=strcat(path1,'\',listing(4).name);
% A=imread(reference);
% B=imread(moving);

B_trans=imtranslate(B,[x(2),x(3)]);
ang=(180/pi)*x(1);
B_reg=imrotate(B_trans,ang);

figure;
imshow(A);
figure;
imshow(B);
figure;
imshow(B_reg);



