function [Gradient] = Compute_Gradient(M)
%输入值是灰度图或者V通道图
 
[dx,dy]=gradient(double(M));
Gradient=sqrt(dx.^2+dy.^2);
Gradient=Gradient.^0.5;
% sobel=fspecial('Sobel');
% hx_sobel=imfilter(M,sobel,'conv','symmetric');
% hy_sobel=imfilter(M',sobel,'conv','symmetric');
% Gradient=hx_sobel+hy_sobel';

end

