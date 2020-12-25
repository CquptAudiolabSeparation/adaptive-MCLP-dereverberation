function [bx_u,bx_z] = derev_MCLP_QRRLS_ADMM(x,wlen,gamma,Lg,T60)
X = stft_multi(x.',wlen,wlen/4);

LBIN = size(X,1);%ffbin? =wlen/2+1
LN = size(X,2);%帧数
LM = size(X,3);%Nmic?信道数

alfa=0.5;%smoothing constant 
Td=0.05;%源信号没有显著自相关系数的时间跨度，常量50ms
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
processed_spec_u=zeros(LBIN, LN, LM);
processed_spec_z=zeros(LBIN, LN, LM);

%******************
% J=zeros(LBIN,LN);
% e=zeros(LBIN,LN);
% e_w=zeros(LBIN,LN);
%******************
rou = 10^3;
iter = 25;
sigema_u = zeros(LBIN, LN, LM);

for k = 1:LBIN 
    k
    G_hat = zeros(LM*Lg, LM);
%     D =  [delta^(0.5) *eye(Lg*LM,Lg*LM) , 1/sqrt(delta) * ones(Lg*LM,2);zeros(1,Lg*LM+LM)];  %2 * eye(Lg*LM+1,Lg*LM+LM);
    D = eye(Lg*LM+1,Lg*LM+LM);
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
            w_2(k, n) = (norm(temp)^2/LM + eps_1)^(p/2 -1); %增加/LM12.15
            for chan = 1:LM 
                sigema_u(k, n, chan) = sqrt(min(sigema_r(k, n, chan) , x_abs_2(chan)));
            end
        end
        %*******************
        if  n<= tau+Lg-1
            if n <= ND
                d_u = reshape(X(k, n, :), LM ,1);
                d_z = reshape(X(k, n, :), LM ,1);
                for chan = 1:LM
                    processed_spec_u(k, n, chan)=d_u(chan);
                    processed_spec_z(k, n, chan)=d_z(chan);
                end
                %***代价函数*初始化*******
%                 e(k,n) = sum(abs(d_u).^2);
%                 e_w(k,n) =  e(k,n);
%                 if n ==1
%                     J(k,n) =  e_w(k,n);
%                 else
%                     J(k,n) = J(k,n-1) + e_w(k,n); 
%                 end
                %**************************   
            else
                x_tau = [zeros(LM*(Lg-n+tau),1);reshape(X(k, 1: n - tau , :), (n-tau)*LM, 1)];
                X_QR = reshape(X(k,n,:),LM,1);
                D(LM*Lg+1 , :) = sqrt(w_2(k, n)) .* [x_tau' , X_QR'];
                R_qr = qrgv_1(D, Lg*LM);%QR分解
                D = sqrt(gamma) * R_qr; 
                R_qr_MLg =  R_qr(1: Lg*LM, 1: Lg*LM);
                d_qr_M = R_qr(1:Lg*LM , Lg*LM+1:Lg*LM+LM);
                for chan = 1:LM
                    G_hat(:,chan) = backsolution(R_qr_MLg , d_qr_M(:,chan));
                end
                              
                u_hat = G_hat' * x_tau;
                           
                Q_hat_i = inv( R_qr_MLg' * R_qr_MLg );
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
                
                %***代价函数********
%                 e(k,n) = sum(abs(d_u).^2);
%                 e_w(k,n) = w_2(k,n)*e(k,n);
%                 J(k,n) = gamma * J(k,n-1)+ e_w(k,n);
                %*******************  
            end            
            %************************  
        else
            %********QRRLS*******
            X_QR = reshape(X(k,n,:),LM,1);
            x_tau = reshape(X(k, n - tau - Lg + 1: n - tau , :), Lg*LM, 1);
%             D(LM*Lg+1 , :) = sqrt(w_2(k, n)) .* [x_tau' , X_QR(1),X_QR(2)];
            D(LM*Lg+1 , :) = sqrt(w_2(k, n)) .* [x_tau' , X_QR'];

            R_qr = qrgv_1(D, Lg*LM);%QR分解
            R_qr_MLg =  R_qr(1: Lg*LM, 1: Lg*LM);
            d_qr_M = R_qr(1:Lg*LM , Lg*LM+1:Lg*LM+LM);
            for chan = 1:LM
                G_hat(:,chan) = backsolution(R_qr_MLg , d_qr_M(:,chan));
            end
                        
            u_hat = G_hat' * x_tau;
            D = sqrt(gamma) * R_qr; 
            
            Q_hat_i = inv( R_qr_MLg' * R_qr_MLg );
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
            
            %********************
        
            %***存储中间变量*****
%             ye=n+(k-1)*LN;
%             G(:,:,ye) = G_hat;
            %*******************
        
            %***代价函数********
%             e(k,n) = sum(abs(d_u).^2);
%             e_w(k,n) = w_2(k,n)*e(k,n);
%             J(k,n) = gamma*J(k,n-1)+ e_w(k,n);
            %*******************      
        end
    end
end
bx_u = istft_multi(processed_spec_u,size(x ,1))';
bx_z = istft_multi(processed_spec_z,size(x ,1))';
end

%save mic_p=0_Lg=20_r=0.965_QR_AMCLP.mat;


