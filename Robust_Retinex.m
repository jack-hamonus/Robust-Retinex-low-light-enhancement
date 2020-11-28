function [R,L,N] = Robust_Retinex(I,G)
%I是要生成梯度图的原图对应的hsv模型下的v通道
%   Micro=µ ，其初始值为1
%    Rho=ρ 其初始值为1.5
%    Beta=β 其值为 0.05
%    Omega=ω 其值为0.01
%    Delta=δ 其值为1
Micro=1;
Rho=1.5;
Beta=0.05;
Omega=0.01;
Delta=1;

[rows,columns]=size(I);

%初始化
L_pre=I;
R_pre=zeros(rows,columns);

T_pre=zeros(rows,columns);
Z_pre=T_pre;
N_pre=T_pre;
     
i=I(:);
g=G(:);

n=N_pre(:);
t=T_pre(:);
z=Z_pre(:);
l_pre=L_pre(:);
l_next=l_pre;
r_pre=R_pre(:);
r_next=r_pre;




%l_diagonal_matrix represents a diagonal matrix with l as its entries.

Length_l=size(l_pre);
Length_l=Length_l(1);




%表示D,这里D的定义仍在存在问题，需要商议
%  Ix = diff(L_pre,1,2); Ix = padarray(Ix, [0 1], 'post');
%  Iy = diff(L_pre,1,1); Iy = padarray(Iy, [1 0], 'post');
% Ix_vector=Ix(:);
% Iy_vector=Iy(:);

% ux = padarray(Ix_vector, rows, 'pre'); ux = ux(1:end-rows);
%     uy = padarray(  Iy_vector, 1, 'pre'); uy = uy(1:end-1);
%     D = Ix_vector+ux+Iy_vector+uy;
%     T = spdiags([-Ix_vector, -Iy_vector],[-rows,-1],Length_l,Length_l);
%     MN = T + T' + spdiags(D, 0, Length_l, Length_l); %MN代表论文中D

%求D新方法 2020、11、20
% MN=zeros(1,Length_l);
% MN(1)=1;
% MN(2)=-1;
% MN=toeplitz(MN);


e=ones(Length_l,1);
MN=spdiags([-e,2*e,-e],[-1,1],Length_l,Length_l);
%K表示收敛算法迭代的次数





for k=1:20
    disp('Iterator: ')
    disp(k);
    
% R sub-problem
  r_pre=r_next;
  l_diagonal_matrix = spdiags(l_next,0,Length_l,Length_l);
  
%   r_next=inv((l_diagonal_matrix*l_diagonal_matrix')+Omega*(D*D'))*(l_diagonal_matrix*(i-n)+Omega*D'*g);
 DEN=l_diagonal_matrix.^2+Omega*MN*MN';
 DEN=sparse(DEN);
 NUM=l_diagonal_matrix*(i-n)+Omega*MN'*g; 
%  A = ichol(DEN,struct('michol','on'));
%  [r_next,~] = pcg(DEN, NUM(:), 0.01, 40,A,A'); 
%  [L,U] = ilu(DEN,struct('type','ilutp','droptol',0.01));
%             [r_next,~] = bicg(DEN,NUM(:), 0.01, 40, L, U);
[r_next,~] = minres(DEN,NUM(:), 0.01, 40);
% r_next=inv(DEN)*NUM;
 R_next=reshape(r_next,rows,columns);   
 disp('R over');
 
%L sub-problem
l_pre=l_next;
r_diagonal_metrix=spdiags(r_next,0,Length_l,Length_l);
DEN=2*r_diagonal_metrix.^2+Micro*MN*MN'; 
DEN=sparse(DEN);
NUM=2*r_diagonal_metrix*(i-n)+Micro*MN'*(t-(z/Micro));
% A = ichol(DEN,struct('michol','on'));
% [l_next,~] = pcg(DEN, NUM(:), 0.01, 40,A,A'); 


% [L,U] = ilu(DEN,struct('type','ilutp','droptol',0.01));
% [l_next,~] = bicg(DEN,NUM(:), 0.01, 40, L, U);

[L_next,~] = minres(DEN,NUM(:), 0.01, 40);

% l_next=inv(DEN)*NUM;
  L_next=reshape(l_next,rows,columns);
disp('L over');
% N sub-problem
N_pre=(1-(R_next.*L_next))/(1+ Delta);
n=N_pre(:);


%T sub-problem
data=Compute_Gradient(L_next)+Z_pre/Micro;
T=sign(data).*max(abs(data)-0.01,0);
t=T(:);
%Update Z aND Micro
Z_pre=Z_pre+Micro*(Compute_Gradient(L_next)-T);
z=Z_pre(:);
Micro=Micro*Rho;
if(norm(r_next-r_pre,'fro')<0.01)
    break;
end






end
R=R_next;
L=L_next;
N=N_pre;

