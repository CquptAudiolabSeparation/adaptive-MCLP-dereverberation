function [bx,e,e_w,J,parm] = derev_MCLP_RLS(x,wlen,gamma,Lg,T60)
X = stft_multi(x.',wlen,wlen/4);

LBIN = size(X,1);%ffbin? =wlen/2+1
LN = size(X,2);%帧数
LM = size(X,3);%Nmic?信道数

alfa=0.5;%smoothing constant 
Td=0.05;%源信号没有显著自相关系数的时间跨度，常量50ms
% T60=0.7;%T60混响时间
deta=3*log(10)/(T60/Td);%decay constant
p=0.5;    % p norm
% ita=10^-8;% regularization parameter
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
J = zeros(LBIN,LN);
e = zeros(LBIN,LN);
e_w = zeros(LBIN,LN);

parm.G = cell(LBIN,LN);
%Gd(:,:,ye) = reshape(G_downhat(:,:,I+1),LM*Lg,LM); 
parm.Q_i = cell(LBIN,LN);
parm.K = cell(LBIN,LN);
%******************

for k = 1:LBIN
    k
    for n=1:  tau + Lg - 1
        d = reshape(X(k, n, :), LM ,1);
        for chan = 1:LM
            processed_spec(k, n, chan)=d(chan);
        end
    end
    %Algorithm 2:PSD Estimation
    for n = 1 + ND:LN
        x_abs_2 = abs(X(k, n, :)).^2;
        sigema_x(k, n, :) = alfa*sigema_x(k, n-1, :) + alfa_1*x_abs_2;
        sigema_r(k, n, :) = coff1*sigema_x(k, n - ND, :);
        sigema_d(k, n, :) = alfa*sigema_d(k, n-1, :) + alfa_1* max((x_abs_2 - sigema_r(k, n, :)), 0);
        temp = reshape(sigema_d(k, n, :), LM,1);
%         w_2(k, n) = (norm(temp)^2 + eps_1)^(p/2 -1); 
        w_2(k, n) = (norm(temp)^2/LM + eps_1)^(p/2 -1); %增加/LM12.15
    end
    
    Q_hat_i = 10^(-4) * eye(Lg*LM);
    Q_hat_i_old = Q_hat_i;

    eye_1 = eye(Lg * LM);
    
    %***代价函数*初始化**
    e(k,1) = sum(abs(X(k, 1, :)).^2);
    e_w(k,1) = w_2(k,1) * e(k,1);
    J(k,1) = e_w(k,1);  
   
    for n= 2 : tau + Lg - 1
    e(k,n) = sum(abs(X(k, n, :)).^2);
    e_w(k,n) = w_2(k,n) * e(k,n);
    J(k,n) =  J(k,n-1) + e_w(k,n);
    end
    %*****************
    
    G_hat_old =  zeros(LM*Lg, LM);
    
    for n = 1 + tau + Lg - 1 :LN
        x_tau = reshape(X(k, n - tau :-1: n - tau - Lg + 1, :), Lg*LM, 1);
        k_hat = (Q_hat_i_old * x_tau)/((gamma / w_2(k, n)) + x_tau' * Q_hat_i_old * x_tau);%dengyu0
        G_hat = G_hat_old +  k_hat * (reshape(X(k, n, :), LM ,1) - G_hat_old' * x_tau)';%bugenxin
        Q_hat_i = (1/ gamma) * (eye_1 - k_hat * x_tau')*Q_hat_i_old;
        
        % update of weights
%         sumG = sum(sum(abs(G_hat).^2));
%         if sumG >Beta
%             G_hat = sqrt(Beta ./ (sumG * ones(Lg*LM,LM))) .* G_hat ;
%         end
        
        u = G_hat' * x_tau;
        
        G_hat_old = G_hat;
        Q_hat_i_old = Q_hat_i;
        d = reshape(X(k, n, :), LM ,1) - u;
        for chan = 1:LM
            processed_spec(k, n, chan)=d(chan);
        end
        
        %***存储中间变量*****
%         ye=n+(k-1)*LN;
%         parm.G{k,n} = G_hat;
%          %Gd(:,:,ye) = reshape(G_downhat(:,:,I+1),LM*Lg,LM); 
%         parm.Q_i{k,n} = Q_hat_i;
%         parm.K{k,n} = k_hat;
        
        %*******************
        
        %***代价函数********
        e(k,n) = sum(abs(d).^2);
        e_w(k,n) = w_2(k,n)*e(k,n);
        J(k,n) = gamma*J(k,n-1)+ e_w(k,n);
        %*******************
    end
end
bx = istft_multi(processed_spec,size(x ,1))';
% bx = bx./repmat(max(abs(bx)) , size(bx,1), 1);
end
% for k = 4:4:LBIN
%     a= figure;plot(1:LN, w_2(k-1, :), 1:LN, w_2(k-2, :), 1:LN, w_2(k-3, :),1:LN, w_2(k, :));
%     close(a);
% end

% save mic_p=0_Lg=20_r=0.965_AMCLP.mat;


