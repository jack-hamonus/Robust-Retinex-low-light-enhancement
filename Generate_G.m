function [G] = Generate_G(v)
%gimg_gray是要生成梯度图的原图对应的hsv模型下的v通道
%   该论文中λ=10，d=10,eps=0.01
                                                                               
[dx,dy]=gradient(v);
G=sqrt(dx.^2+(dy.^2));
% v=double(v);
% sobel=fspecial('Sobel');
% hx_sobel=imfilter(v,sobel,'conv','symmetric');
% hy_sobel=imfilter(v',sobel,'conv','symmetric');
% G=hx_sobel+hy_sobel';
G(G<eps)=0;
K=1+10*exp(-abs(G)/10);
G=K.*G;
end

