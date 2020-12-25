function x=backsolution(A,b)
% 求解 上三角类型 的线性方程组,
% 使用回代法
% 于 2015/11/10
% yanghsj@ahut.edu.cn
n=length(b);   %右端项
x=zeros(n,1);  %解为列向量
x(n)=b(n)/A(n,n);
for i=n-1:-1:1
    x(i)=(b(i)-A(i,i+1:n)*x(i+1:n))/A(i,i);
end