function [Gradient] = Compute_Gradient(M)
%����ֵ�ǻҶ�ͼ����Vͨ��ͼ
 
[dx,dy]=gradient(double(M));
Gradient=sqrt(dx.^2+dy.^2);
Gradient=Gradient.^0.5;
% sobel=fspecial('Sobel');
% hx_sobel=imfilter(M,sobel,'conv','symmetric');
% hy_sobel=imfilter(M',sobel,'conv','symmetric');
% Gradient=hx_sobel+hy_sobel';

end

