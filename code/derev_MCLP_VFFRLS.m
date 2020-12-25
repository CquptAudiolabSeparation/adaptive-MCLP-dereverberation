function [bx gam] = derev_MCLP_VFFRLS(x,wlen,gamma_lower,gamma_upper,Lg,T60)
X = stft_multi(x.',wlen,wlen/4);
LBIN = size(X,1);%ffbin? =wlen/2+1
LN = size(X,2);%帧数
LM = size(X,3);%Nmic?信道数

alfa=0.5;%smoothing constant 
Td=0.05;%源信号没有显著自相关系数的时间跨度，常量50ms
deta=3*log(10)/(T60/Td);%decay constant
p=0.5;    % p norm

coff1 = exp(-2*deta);
% ND = floor(Td/0.016);
iota = 0.99;
T = Lg;
M = 20;

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
% J=zeros(LBIN,LN);
% e=zeros(LBIN,LN);
% e_w=zeros(LBIN,LN);

G = zeros(Lg*LM,LM,LN*LBIN); %G_hat
gam = zeros(LBIN,LN);

Beta = 5 * LM;

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
%        w_2(k, n) = (norm(temp)^2 + eps_1)^(p/2 -1); 
        w_2(k, n) = (norm(temp)^2/LM + eps_1)^(p/2 -1); %增加/LM12.15
    end
    
    Q_hat_i = 10^(-8) * eye(Lg*LM);
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
%     J(k,n) = gamma_upper * J(k,n-1) + e_w(k,n);
%     end
    %*******************
    
    G_hat_old = 1e-10* ones(LM*Lg, LM);
    
    G_VFF_old = 1e-10* ones(LM*Lg, LM);
    C_VFF_old = 1e-10* ones(LM*Lg, LM);
    grad_C = zeros(LN, 1);  
    grad_C_ave = zeros(LN, 1); 
    grad_N = zeros(LN, 1);
    
    for n = 1 + tau + Lg - 1 :LN
        %********VFF*********
        G_VFF = iota * G_VFF_old + (1 - iota) * G_hat_old; 
        C_VFF = G_hat_old - G_VFF_old;
        grad_C(n,:) = abs(norm(C_VFF,1) - norm(C_VFF_old,1));
            
        if n>(tau + Lg + T - 2)
            grad_C_ave(n,:) = mean(grad_C(n-T+1:n,:));
        else
            grad_C_ave(n,:) = mean(grad_C(tau + Lg : n,:));
        end
        
        if n == tau + Lg
            grad_N(n,:) = 1;
        elseif n <= tau + Lg + M
            grad_N(n,:) = (grad_C_ave(n,:) - min(grad_C_ave(tau + Lg :n,:))) / (max(grad_C_ave(tau + Lg :n,:)) - min(grad_C_ave(tau + Lg :n,:)));
        else
            grad_N(n,:) = min(abs(grad_C_ave(n,:) ./ grad_C_ave(n - M,:)),1);
%             grad_N(n,:) = (grad_C_ave(n,:) - min(grad_C_ave(n - M + 1 :n,:))) / (max(grad_C_ave(n - M + 1 :n,:)) - min(grad_C_ave(n - M + 1 :n,:)));

        end
%         if n > (tau + Lg + M - 2)
%             grad_N(n,:) = min(abs(grad_C_ave(n,:) ./ grad_C_ave(tau + Lg + M - 2,:)),1);
%         else
%             grad_N(n,:) = 1;
%         end
        gamma = gamma_lower + (1 - grad_N(n,:)) * (gamma_upper - gamma_lower);
        %********************
        
        x_tau = reshape(X(k, n - tau :-1: n - tau - Lg + 1, :), Lg*LM, 1);
        k_hat = (Q_hat_i_old * x_tau)/((gamma / w_2(k, n)) + x_tau' * Q_hat_i_old * x_tau);%dengyu0
        G_hat = G_hat_old +  k_hat * (reshape(X(k, n, :), LM ,1) - G_hat_old' * x_tau)';%bugenxin
        Q_hat_i = (1/ gamma) * (eye_1 - k_hat * x_tau')*Q_hat_i_old;
        
        % update of weights
        sumG = sum(sum(abs(G_hat).^2));
        if sumG >Beta
            G_hat = sqrt(Beta ./ (sumG * ones(Lg*LM,LM))) .* G_hat ;
        end
        
        u = G_hat' * x_tau;
        
        G_hat_old = G_hat;
        Q_hat_i_old = Q_hat_i;
        d = reshape(X(k, n, :), LM ,1) - u;
        for chan = 1:LM
            processed_spec(k, n, chan)=d(chan);
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
        
        %******VFF*********
        gam(k,n) = gamma;
        G_VFF_old = G_VFF;
        C_VFF_old = C_VFF;
        %******************
        
    end
end
bx = istft_multi(processed_spec,size(x ,1))';
end
% for k = 4:4:LBIN
%     a= figure;plot(1:LN, w_2(k-1, :), 1:LN, w_2(k-2, :), 1:LN, w_2(k-3, :),1:LN, w_2(k, :));
%     close(a);
% end

% save mic_p=0_Lg=20_0.9_0.99_VFF_AMCLP.mat;


