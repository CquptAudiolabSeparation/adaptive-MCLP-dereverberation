function R=qrgv_1(A,num)
% ����Givens�任��������A�ֽ�ΪA=QR������QΪ��������RΪ��������
%
% ����˵��
% A����Ҫ����QR�ֽ�ķ���
% R���ֽ�õ�����������
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




