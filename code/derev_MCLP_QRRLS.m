function [bx,e,e_w,J,parm] = derev_MCLP_QRRLS(x,wlen,gamma,Lg,T60)
X = stft_multi(x.',wlen,wlen/4);

LBIN = size(X,1);%ffbin? =wlen/2+1
LN = size(X,2);%֡��
LM = size(X,3);%Nmic?�ŵ���

alfa=0.5;%smoothing constant 
Td=0.05;%Դ�ź�û�����������ϵ����ʱ���ȣ�����50ms
deta=3*log(10)/(T60/Td);%decay constant
p=1/2;    % p norm
% ita=10^-8;% regularization parameter
% gamma=0.75; %forgetting parameter

coff1 = exp(-2*deta);
% ND = floor(Td/0.016);

ND = 4;
tau = 2;

sigema_d = zeros(LBIN, LN, LM);
sigema_x = zeros(LBIN, LN, LM);
sigema_r = zeros(LBIN, LN, LM);
w_2 = zeros(LBIN, LN);

alfa_1 = (1-alfa);
eps_1 = eps*10;
processed_spec=zeros(LBIN, LN, LM);
%******************
J=zeros(LBIN,LN);
e=zeros(LBIN,LN);
e_w=zeros(LBIN,LN);
%******************
Beta = 5 * LM;
parm.G = cell(LBIN,LN);
parm.R_qr_MLg = cell(LBIN,LN);
parm.d_qr_M = cell(LBIN,LN);

for k = 1:LBIN 
    k
    G_hat = zeros(LM*Lg, LM);
%     D =  [delta^(0.5) *eye(Lg*LM,Lg*LM) , 1/sqrt(delta) * ones(Lg*LM,2);zeros(1,Lg*LM+LM)];  %2 * eye(Lg*LM+1,Lg*LM+LM);
    D = 10^2  * eye(Lg*LM+1,Lg*LM+LM);
%     D = [triu(ones(LM*Lg)),ones(Lg*LM,LM);zeros(1,Lg*LM+LM)];
%     D = zeros(Lg*LM+1,Lg*LM+2);

    for n = 1:LN  
        if n >= 1+ND
            %Algorithm 2:PSD Estimation
            x_abs_2 = abs(X(k, n, :)).^2;
            sigema_x(k, n, :) = alfa*sigema_x(k, n-1, :) + alfa_1*x_abs_2;
            sigema_r(k, n, :) = coff1*sigema_x(k, n - ND, :);
            sigema_d(k, n, :) = alfa*sigema_d(k, n-1, :) + alfa_1* max((x_abs_2 - sigema_r(k, n, :)), 0);
            temp = reshape(sigema_d(k, n, :), LM,1);
            w_2(k, n) = (norm(temp)^2/LM + eps_1)^(p/2 -1); %����/LM12.15
        end
        %*******************
        if  n<= tau+Lg-1
%             if n <= ND
                d = reshape(X(k, n, :), LM ,1);
                for chan = 1:LM
                    processed_spec(k, n, chan)=d(chan);
                end
                %***���ۺ���*��ʼ��*******
                e(k,n) = sum(abs(d).^2);
                e_w(k,n) =  w_2(k,n) * e(k,n);
                if n ==1
                    J(k,n) =  e_w(k,n);
                else
                    J(k,n) = J(k,n-1) + e_w(k,n); 
                end
                %**************************   
        else
            %********QRRLS*******
            X_QR = reshape(X(k,n,:),LM,1);
            x_tau = reshape(X(k, n - tau - Lg + 1: n - tau , :), Lg*LM, 1);
%             D(LM*Lg+1 , :) = sqrt(w_2(k, n)) .* [x_tau' , X_QR(1),X_QR(2)];
            D(LM*Lg+1 , :) = sqrt(w_2(k, n)) .* [x_tau' , X_QR'];

            R_qr = qrgv_1(D, Lg*LM);%QR�ֽ�
            R_qr_MLg =  R_qr(1: Lg*LM, 1: Lg*LM);
            d_qr_M = R_qr(1:Lg*LM , Lg*LM+1:Lg*LM+LM);
            for chan = 1:LM
                G_hat(:,chan) = backsolution(R_qr_MLg , d_qr_M(:,chan));
            end
            
            % update of weights
%             sumG = sum(sum(abs(G_hat).^2));
%             if sumG >Beta
%                 G_hat = sqrt(Beta ./ (sumG * ones(Lg*LM,LM))) .* G_hat ;
%             end
            
            u = G_hat' * x_tau;
            d = reshape(X(k, n, :), LM ,1) - u;
            for chan = 1:LM
                processed_spec(k, n, chan)=d(chan);
            end
            
            D = sqrt(gamma) * R_qr; 
            %********************
        
            %***�洢�м����*****
%             ye=n+(k-1)*LN;
%             parm.G{k,n} = G_hat;
%             parm.R_qr_MLg{k,n} = R_qr_MLg;
%             parm.d_qr_M{k,n} = d_qr_M;
            
            %*******************
        
            %***���ۺ���********
            e(k,n) = sum(abs(d).^2);
            e_w(k,n) = w_2(k,n)*e(k,n);
            J(k,n) = gamma*J(k,n-1)+ e_w(k,n);
            %*******************      
        end
    end
end
bx = istft_multi(processed_spec,size(x ,1))';
% bx = bx./repmat(max(abs(bx)) , size(bx,1), 1);
end

%save mic_p=0_Lg=20_r=0.965_QR_AMCLP.mat;


