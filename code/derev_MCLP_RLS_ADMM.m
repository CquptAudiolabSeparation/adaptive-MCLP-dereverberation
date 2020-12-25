function [bx_u,bx_z] = derev_MCLP_RLS_ADMM(x,wlen,gamma,Lg,T60)
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
processed_spec_u=zeros(LBIN, LN, LM);
processed_spec_z=zeros(LBIN, LN, LM);
%******************
% J = zeros(LBIN,LN);
% e = zeros(LBIN,LN);
% e_w = zeros(LBIN,LN);
%******************
rou = 10^3;
iter = 25;
sigema_u = zeros(LBIN, LN, LM);

for k = 1:LBIN
    k
    for n=1:  tau + Lg - 1
        d_u = reshape(X(k, n, :), LM ,1);
        d_z = reshape(X(k, n, :), LM ,1);
        for chan = 1:LM
            processed_spec_u(k, n, chan)=d_u(chan);
            processed_spec_z(k, n, chan)=d_z(chan);
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
        %ADMM
        for chan = 1:LM 
            sigema_u(k, n, chan) = sqrt(min(sigema_r(k, n, chan) , x_abs_2(chan)));
        end
    end
    
    Q_hat_i = 10^(-4) * eye(Lg*LM);
    Q_hat_i_old = Q_hat_i;

    eye_1 = eye(Lg * LM);
    
    %***代价函数*初始化**
%     e(k,1) = sum(abs(X(k, 1, :)).^2);
%     e_w(k,1) = w_2(k,1) * e(k,1);
%     J(k,1) = e_w(k,1);  
%    
%     for n= 2 : tau + Lg - 1
%     e(k,n) = sum(abs(X(k, n, :)).^2);
%     e_w(k,n) = w_2(k,n) * e(k,n);
%     J(k,n) = gamma * J(k,n-1) + e_w(k,n);
%     end
    %*****************
    
    G_hat_old =  zeros(LM*Lg, LM);
    
    for n = 1 + tau + Lg - 1 :LN
        %algorithm 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x_tau = reshape(X(k, n - tau :-1: n - tau - Lg + 1, :), Lg*LM, 1);
        k_hat = (Q_hat_i_old * x_tau)/((gamma / w_2(k, n)) + x_tau' * Q_hat_i_old * x_tau);%dengyu0
        G_hat = G_hat_old +  k_hat * (reshape(X(k, n, :), LM ,1) - G_hat_old' * x_tau)';%bugenxin
        Q_hat_i = (1/ gamma) * (eye_1 - k_hat * x_tau')*Q_hat_i_old;
                
        u_hat = G_hat' * x_tau;
        G_hat_old = G_hat;
        Q_hat_i_old = Q_hat_i;
        
        %algorithm 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        z_conhat = zeros(LM,1);lamda = zeros(LM,1);
        k_conhat = (Q_hat_i * x_tau)/(2 / rou + x_tau' * Q_hat_i * x_tau);
        for it = 1:iter
            G_conhat = G_hat + k_conhat * (z_conhat + lamda - u_hat)';
            u_conhat = G_conhat' * x_tau;
            for chan = 1:LM
                z_conhat(chan,1) = min(sigema_u(k, n, chan)/abs(u_conhat(chan) - lamda(chan)) , 1) * (u_conhat(chan) - lamda(chan));
            end
            lamda = lamda + z_conhat - u_conhat;
        end
        d_u = reshape(X(k, n, :), LM ,1) - u_conhat;
        d_z = reshape(X(k, n, :), LM ,1) - z_conhat;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for chan = 1:LM
            processed_spec_u(k, n, chan)=d_u(chan);
            processed_spec_z(k, n, chan)=d_z(chan);
        end
        
        %***存储中间变量*****
%         ye=n+(k-1)*LN;
%         G(:,:,ye) = G_hat;
%          %Gd(:,:,ye) = reshape(G_downhat(:,:,I+1),LM*Lg,LM); 
%         Q_i(1,:,ye) = diag(Q_hat_i);
%         K(:,:,ye) = k_hat;
         %Kd(:,:,ye) = k_downhat;
        %*******************
        
        %***代价函数********
%         e(k,n) = sum(abs(d).^2);
%         e_w(k,n) = w_2(k,n)*e(k,n);
%         J(k,n) = gamma*J(k,n-1)+ e_w(k,n);
        %*******************
    end
end
bx_u = istft_multi(processed_spec_u,size(x ,1))';
% bx_u = bx_u./repmat(max(abs(bx_u)) , size(bx_u,1), 1);

bx_z = istft_multi(processed_spec_z,size(x ,1))';
% bx_z = bx_z./repmat(max(abs(bx_z)) , size(bx_z,1), 1);

end
% for k = 4:4:LBIN
%     a= figure;plot(1:LN, w_2(k-1, :), 1:LN, w_2(k-2, :), 1:LN, w_2(k-3, :),1:LN, w_2(k, :));
%     close(a);
% end

% save mic_p=0_Lg=20_r=0.965_AMCLP.mat;


