%This function computes mutual information (MI) between two images (reference image and moving image) 
% and returns -MI. Input parameter x determines the angle of rotation, and row-wise
% and column-wise translation. The code is fully vectorized and does not
% use any loop

function[MI_n2]=Mutual_Information_Computation(x)

%Hard coded since there are only two images
A=imread('Data\Face1.jpg'); %Reference Image
B=imread('Data\Face1_Rotate_15_X250_Y-24.jpg'); %Moving Image

%Fetching data from location

% path1=strcat(pwd,'\Data');
% listing=dir(path1);
% reference=strcat(path1,'\',listing(3).name);
% moving=strcat(path1,'\',listing(4).name);
% A=imread(reference);
% B=imread(moving);

%Translation
B=imtranslate(B,[x(2),x(3)]);

A=double(A); %Converting unsigned integer values to double precision values
B=double(B);

%Coordinate Matrix Formation
s=size(B); %Image B contains s(1) rows and s(2) columns

%Center coordinates of B are tx and ty
tx=round(s(1)/2); 
ty=round(s(2)/2);

%Coordinate matrices are formed for rotation
[corx,cory]=ndgrid(1:s(1),1:s(2));
corx1=reshape(corx,[1,s(1)*s(2)]);
cory1=reshape(cory,[1,s(1)*s(2)]);
lr=ones(1,s(1)*s(2));
cor_mat=[corx1;cory1;lr];

%Translation Matrix Formation
tm=[1,0,-tx;0,1,-ty;0,0,1];

%First Transformation = Translation
TC=tm*cor_mat;

%Rotation Matrix Formation. 'x' is rotation angle
cosx=cos(x(1));
sinx=sin(x(1));
rot_mat=[cosx,-sinx,0;sinx,cosx,0;0,0,1];

%Second Transformation = Rotation
rot_cor=rot_mat*TC;

%Reshaping transformed image matrix
rot_mov_row=reshape(rot_cor(1,:),[s(1) s(2)]);
rot_mov_col=reshape(rot_cor(2,:),[s(1) s(2)]);

%Center coordinates of A
s2=size(A);
tx2=round(s2(1)/2);
ty2=round(s2(2)/2);

%Coordinate matrices for A
corx2=repmat((1:s2(1))',1,s2(2));
cory2=repmat((1:s2(2)),s2(1),1);
corx22=reshape(corx2,[1,s2(1)*s2(2)]);
cory22=reshape(cory2,[1,s2(1)*s2(2)]);
lr2=ones(1,s2(1)*s2(2));
cor_mat2=[corx22;cory22;lr2];

%Translation of center of A to (0,0)
tm2=[1,0,-tx2;0,1,-ty2;0,0,1];
TC2=tm2*cor_mat2;

trns_trf_row=reshape(TC2(1,:),[s2(1) s2(2)]);
trns_trf_col=reshape(TC2(2,:),[s2(1) s2(2)]);

%Finding the extreme corner points of A, after translation
min_ref_row=min(trns_trf_row(:));
max_ref_row=max(trns_trf_row(:));
min_ref_col=min(trns_trf_col(:));
max_ref_col=max(trns_trf_col(:));

%Vectorization starts

%Finding floor and ceiling values of rotated coordinates of B
mov_x_flV=floor(rot_mov_row);
mov_x_clV=ceil(rot_mov_row);
mov_y_flV=floor(rot_mov_col);
mov_y_clV=ceil(rot_mov_col);

%Finding the coordinates which are whole numbers and making ceiling=floor+1
%for those coordinates

diff_row=mov_x_flV-mov_x_clV;
XX=(find(diff_row==0));
[row,col]=ind2sub(size(diff_row),XX);
mov_x_clV(row,col)=mov_x_flV(row,col)+1;

diff_col=mov_y_flV-mov_y_clV;
YY=find(diff_col==0);
[row,col]=ind2sub(size(diff_col),YY);
mov_y_clV(row,col)=mov_y_flV(row,col)+1;

%Finding out distance of rotated points from nearest ceiling and floor values in order
%to compute weights

difx1V=rot_mov_row-mov_x_flV;
difx2V=1-difx1V;
dify1V=rot_mov_col-mov_y_flV;
dify2V=1-dify1V;

%Vectorized weight computation

w11V=difx2V.*dify2V;
w12V=difx2V-w11V;
w21V=dify2V-w11V;
w22V=difx1V-w21V;

%Vectorized Joint Histogram

%Finding out the pixels in the rotated image, which lie within the
%boundaries of the reference image (A). In other words, finding out the
%overlapped region

Pos11=find(mov_x_flV>=min_ref_row & mov_x_flV<=max_ref_row & mov_y_flV>=min_ref_col & mov_y_flV<=max_ref_col);
Pos12=find(mov_x_flV>=min_ref_row & mov_x_flV<=max_ref_row & mov_y_clV>=min_ref_col & mov_y_clV<=max_ref_col);
Pos21=find(mov_x_clV>=min_ref_row & mov_x_clV<=max_ref_row & mov_y_flV>=min_ref_col & mov_y_flV<=max_ref_col);
Pos22=find(mov_x_clV>=min_ref_row & mov_x_clV<=max_ref_row & mov_y_clV>=min_ref_col & mov_y_clV<=max_ref_col);

%Finding out values of mov_x_flV, mov_x_clV, mov_y_flV and mov_y_clV for
%indices Pos11,Pos12, Pos21 and Pos22 and adding center coordinates

X_Pos=[mov_x_flV(Pos11);mov_x_flV(Pos12);mov_x_clV(Pos21);mov_x_clV(Pos22)];
X_Pos=X_Pos+tx2;
Y_Pos=[mov_y_flV(Pos11);mov_y_clV(Pos12);mov_y_flV(Pos21);mov_y_clV(Pos22)];
Y_Pos=Y_Pos+ty2;

%Converting X_Pos and Y_Pos into linear indices
lin_Ind=sub2ind(size(A),X_Pos,Y_Pos);

%Intensities of image A, at indices lin_Ind
R=A(lin_Ind);

%Finding out intensities of moving image(B), at positions Pos11, Pos12,
%Pos21 and Pos22

Pos=[Pos11;Pos12;Pos21;Pos22];
M=B(Pos);

%Finding out weights at positions Pos11, Pos12, Pos21 and Pos22

w11P=w11V(Pos11);
w12P=w12V(Pos12);
w21P=w21V(Pos21);
w22P=w22V(Pos22);

%Formation of joint histogram with values W, to be placed at positions T
W=[w11P;w12P;w21P;w22P];
T=[M+1,R+1];
hist_J= accumarray(T,W,[],@(x) sum(x,'native'));

%%%%%%%%%%%%%%%%%%%%%%%% PROPOSED FORMULA BEGINS %%%%%%%%%%%%%%%%%%%%%%%%%%

% % %Entropy and Joint Entropy Computation with the proposed formulae

% S_hist_J=size(hist_J);
% diff_r=256-S_hist_J(1);
% diff_c=256-S_hist_J(2);
% hist_J=padarray(hist_J,[diff_r,diff_c],'post');
% 
% %Normalization of histogram
% SS=sum(hist_J(:));
% hist_J=hist_J/SS;
%
% Row wise and column wise summation of joint histogram gives individual histograms of A and B 
% hm=sum(hist_J,1);
% hm=hm+eps;
% hr=sum(hist_J,2);
% hr=hr+eps;
% hr=hr';
%
% ent_R1=(-1)*sum(log2((hr(1:2:end-1)+hr(2:2:end))/2).*(hr(1:2:end-1)+hr(2:2:end)));
% ent_M1=(-1)*sum(log2((hm(1:2:end-1)+hm(2:2:end))/2).*(hm(1:2:end-1)+hm(2:2:end)));
% A=hist_J(1:2:end-1,1:2:end-1);
% B=hist_J(1:2:end-1,2:2:end);
% C=hist_J(2:2:end,2:2:end);
% D=hist_J(2:2:end,1:2:end-1);
% X=A+B+C+D;
% Y=X.*log2(X/4+eps);
% ent_J1=(-1)*sum(sum(Y));

%%%%%%%%%%%%%%%%%%%%%%%%%%% PROPOSED FORMULA ENDS %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% STANDARD FORMULA BEGINS %%%%%%%%%%%%%%%%%%%%%%%

%Entropy and Joint Entropy Computation with the original formulae

%Normalization of histogram
SS=sum(hist_J(:));
hist_J=hist_J/SS;

%Row wise and column wise summation of joint histogram gives individual histograms of A and B 
hm=sum(hist_J,1);
hm=hm+eps;
hr=sum(hist_J,2);
hr=hr+eps;
hr=hr';

ent_R1=(-1)*sum(log2((hr(1:end))).*(hr(1:end)));
ent_M1=(-1)*sum(log2((hm(1:end))).*(hm(1:end)));

hist_J=hist_J+eps;
ent_J1=(-1)*sum(sum(log2(hist_J(1:end,1:end)).*(hist_J(1:end,1:end))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STANDARD FORMULA ENDS %%%%%%%%%%%%%%%%%%%%

%Mutual Information Computation 
MI=ent_R1+ent_M1-ent_J1;     
MI_n2=MI*(-1);










