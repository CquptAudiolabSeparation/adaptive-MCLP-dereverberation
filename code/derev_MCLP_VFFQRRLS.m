function [bx,gam, e ,e_w,J,parm] = derev_MCLP_VFFQRRLS(x,wlen,gamma_lower,gamma_upper,Lg,T60)
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
T = 5;
M = 375;

ND = 4;
tau = 2;
delta = 3;

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

parm.G =cell(LBIN,LN); %G_hat
gam = zeros(LBIN,LN);

Beta = 5 * LM;

%******************
for k = 1:LBIN
    k
    G_hat_old =  zeros(LM*Lg, LM);
    G_hat = zeros(LM*Lg, LM);
    
    G_VFF_old = zeros(LM*Lg, LM);
    C_VFF_old = zeros(LM*Lg, LM);
    grad_C = zeros(LN, 1);  
    grad_C_ave = zeros(LN, 1); 
    grad_N = zeros(LN, 1);
    
%     D =  [delta^(0.5) *eye(Lg*LM,Lg*LM) , 1/sqrt(delta) * ones(Lg*LM,2);zeros(1,Lg*LM+LM)];  %2 * eye(Lg*LM+1,Lg*LM+LM);
    D = 10^2 * eye(Lg*LM+1,Lg*LM+LM);
%     D = [triu(ones(LM*Lg)),ones(Lg*LM,2);zeros(1,Lg*LM+LM)];
    
    for n = 1:LN
        if n >= 1 + ND
           %Algorithm 2:PSD Estimation
            x_abs_2 = abs(X(k, n, :)).^2;
            sigema_x(k, n, :) = alfa*sigema_x(k, n-1, :) + alfa_1*x_abs_2;
            sigema_r(k, n, :) = coff1*sigema_x(k, n - ND, :);
            sigema_d(k, n, :) = alfa*sigema_d(k, n-1, :) + alfa_1* max((x_abs_2 - sigema_r(k, n, :)), 0);
            temp = reshape(sigema_d(k, n, :), LM,1);
            w_2(k, n) = (norm(temp)^2/LM + eps_1)^(p/2 -1); %增加/LM12.15  
        end
        if n <=  tau + Lg - 1
%             if n <= ND
                d = reshape(X(k, n, :), LM ,1);
                for chan = 1:LM
                    processed_spec(k, n, chan)=d(chan);
                end 
                %***代价函数*初始化*******
                e(k,n) = sum(abs(d).^2);
                e_w(k,n) =  w_2(k,n) * e(k,n);
                if n ==1
                    J(k,n) =  e_w(k,n);
                else
                    J(k,n) = J(k,n-1) + e_w(k,n); 
                end
                %**************************    

        else
            %********VFF*********
            G_VFF = iota * G_VFF_old + (1 - iota) * G_hat_old; 
            C_VFF = G_hat_old - G_VFF_old;
            grad_C(n,:) = abs(norm(C_VFF,1) - norm(C_VFF_old,1));
                        
            if n>(tau + Lg + T - 2)
                grad_C_ave(n,:) = mean(grad_C(n-T+1:n,:));
            else
                grad_C_ave(n,:) = mean(grad_C(tau + Lg : n,:));       
            end
           
            if n <= tau + Lg + M
                grad_N(n,:) = 1;
            else
%                 if n <= 750
                    tmp = grad_C_ave(n,:) / max(grad_C_ave(n-35:n,:));
%                 else 
%                     tmp = grad_C_ave(n,:) / max(grad_C_ave(n-200:n,:));
%                 end
%                 tmp = grad_C_ave(n,:) / grad_C0;
                if tmp > 1
                    grad_N(n,:)=1;
                   
                elseif tmp < 0
                    grad_N(n,:) = 0;
                else
                    grad_N(n,:) = tmp;
                end

%                 grad_N(n,:) = min(abs(grad_C_ave(n,:) ./ grad_C_ave(tau + Lg + M,:)),1);
%                 grad_N(n,:) = (grad_C_ave(n,:) - min(grad_C_ave(n - M + 1 :n,:))) / (max(grad_C_ave(n - M + 1 :n,:)) - min(grad_C_ave(n - M + 1 :n,:)));
%                   grad_N(n,:) = 1;
            end
%             if grad_N(n,:)  >0.5
%                 h = 2;
%             else 
%                 h = 0.5;
%             end
            h = 1;
            
            
%             if n > (tau + Lg + M - 2)
% %                 grad_N(n,:) = (grad_C_ave(n,:) - min(grad_C_ave(tau + Lg :n,:))) / (max(grad_C_ave(tau + Lg :n,:)) - min(grad_C_ave(tau + Lg :n,:)));
%                 grad_N(n,:) = min(abs(grad_C_ave(n,:) ./ grad_C_ave(tau + Lg + M - 2,:)),1);
%             else
%                 grad_N(n,:) = 1;
%             end
            gamma = gamma_lower + (1 - grad_N(n,:)).^h * (gamma_upper - gamma_lower);
            %********************
            D = sqrt(gamma) * D;
            %********QRRLS*******
            X_QR = reshape(X(k,n,:),LM,1);
            x_tau = reshape(X(k, n - tau - Lg + 1: n - tau , :), Lg*LM, 1);
            D(LM*Lg+1 , :) = sqrt(w_2(k, n)) .* [x_tau' , X_QR'];
            R_qr = qrgv_1(D, Lg*LM);%QR分解
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
            %********************
            D = R_qr; 
            %***存储中间变量*****
%             ye=n+(k-1)*LN;
            parm.G{k,n} = G_hat;
            %*******************
        
            %***代价函数********
            e(k,n) = sum(abs(d).^2);
            e_w(k,n) = w_2(k,n)*e(k,n);
            J(k,n) = gamma*J(k,n-1)+ e_w(k,n);
            %*******************  
        
            %******VFF*********
            gam(k,n) = gamma;
            G_VFF_old = G_VFF;
            C_VFF_old = C_VFF;
            G_hat_old = G_hat;
           %******************
        end
    end
end
bx = istft_multi(processed_spec,size(x ,1))';
% bx = bx./repmat(max(abs(bx)) , size(bx,1), 1);

end

%save mic_p=0_Lg=20_0.9_0.99_VFF_QR_AMCLP.mat;


