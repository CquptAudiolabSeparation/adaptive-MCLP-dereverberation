function x=backsolution(A,b)
% ��� ���������� �����Է�����,
% ʹ�ûش���
% �� 2015/11/10
% yanghsj@ahut.edu.cn
n=length(b);   %�Ҷ���
x=zeros(n,1);  %��Ϊ������
x(n)=b(n)/A(n,n);
for i=n-1:-1:1
    x(i)=(b(i)-A(i,i+1:n)*x(i+1:n))/A(i,i);
end