function R=qrgv_1(A,num)
% 基于Givens变换，将方阵A分解为A=QR，其中Q为正交矩阵，R为上三角阵
%
% 参数说明
% A：需要进行QR分解的方阵
% R：分解得到的上三角阵
R=A;    

for i=1:num
    givens = eye(num+1);
    r = sqrt(abs(R(i,i))^2+abs(R(num+1,i))^2);
    cost=R(i,i)/r;
    sint=R(num+1,i)/r;
    givens(i,i) = conj(cost);
    givens(i,num+1)=conj(sint);
    givens(num+1,i)=-sint;
    givens(num+1,num+1)=cost;     
    R = givens * R;
end




