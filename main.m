addpath('D:\matlab_Demo\Robust_Retinex_code');
img=imread('03.jpg');
img=imresize(img,0.4,'nearest');

img_hsv=rgb2hsv(img);
I=img_hsv(:,:,3);
G=Generate_G(I);
imshow(G);
[R,L,N]=Robust_Retinex(I,G);
% R=real(R);
% L=real(L);
a=1/2.2;
L=L.^a;
% img_v=L.*R;
img_hsv(:,:,3)=L;
img_result=hsv2rgb(img_hsv);
subplot(1,2,1);
imshow(img);
subplot(1,2,2);
imshow(img_result);

% img_v=imresize(img_v,[160,160],'nearest'); 
% l_pre=img_v(:);
% Length_l=size(l_pre);
% Length_l=Length_l(1);
% l_diagonal_matrix=zeros(Length_l,Length_l);
% for temp=1:Length_l
%     l_diagonal_matrix(temp,temp)=l_pre(temp);
% end
% D=zeros(Length_l,Length_l);
% for temp=1:Length_l-1
%     D(temp,temp)=1;
%     D(temp,temp+1)=-1;
% end