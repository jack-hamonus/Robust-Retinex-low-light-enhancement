img=imread('road.jpg');
img=imresize(img,0.7,'nearest');
img_gray=double(rgb2gray(img));
[rows,columns]=size(img_gray);
[dy,dx]=gradient(img_gray);
img_gradient=sqrt(dx.^2+dy.^2);
subplot(1,2,1);
axis();imshow(img_gradient);


%使用Sobel算子
sobel=fspecial('Sobel');
hx=imfilter(img_gray,sobel,'symmetric');
hy=imfilter(img_gray',sobel,'symmetric');
% imshow(hx);
% imshow(hy);
img_grad_sobel=hx+hy';
subplot(1,2,2);
imshow(img_grad_sobel);
grad_matrix=double(img_grad_sobel)./img_gray;

grad_diag_matrix=spdiags(grad_matrix(:),0,rows*columns,rows*columns);


MN=zeros(1,5);
MN(1)=1;
MN(2)=-1;
MN=spdiags();


