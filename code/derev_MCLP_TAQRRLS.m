function [bx,e,e_w,J,parm] = derev_MCLP_TAQRRLS(x,wlen,gamma,Lg,T60)
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
    delta = 10^2 * ones(Lg*LM,1);
    theta = zeros(Lg*LM,1);
    theta1 = zeros(Lg*LM,1);
    
    alpah = 1;
    
%     D = 10^2  * eye(Lg*LM+1,Lg*LM+LM);

    for n = 1 :LN  
        if n >= 1+ND
            %Algorithm 2:PSD Estimation
            x_abs_2 = abs(X(k, n, :)).^2;
            sigema_x(k, n, :) = alfa*sigema_x(k, n-1, :) + alfa_1*x_abs_2;
            sigema_r(k, n, :) = coff1*sigema_x(k, n - ND, :);
            sigema_d(k, n, :) = alfa*sigema_d(k, n-1, :) + alfa_1* max((x_abs_2 - sigema_r(k, n, :)), 0);
            temp = reshape(sigema_d(k, n, :), LM,1);
            w_2(k, n) = (norm(temp)^2/LM + eps_1)^(p/2 -1); %增加/LM12.15
        end
        %*******************
%         if  n<= tau+Lg-1
% %             if n <= ND
%                 d = reshape(X(k, n, :), LM ,1);
%                 for chan = 1:LM
%                     processed_spec(k, n, chan)=d(chan);
%                 end
%                 %***代价函数*初始化*******
%                 e(k,n) = sum(abs(d).^2);
%                 e_w(k,n) =  w_2(k,n) * e(k,n);
%                 if n ==1
%                     J(k,n) =  e_w(k,n);
%                 else
%                     J(k,n) = J(k,n-1) + e_w(k,n); 
%                 end
%         else
            %********QRRLS*******
            if n >= Lg* LM
            X_QR = reshape(X(k,n,:),LM,1);
            x_tau = reshape(X(k, n - tau - Lg + 1: n - tau , :), Lg*LM, 1);
            
            b = - theta;
            b1 = - theta1;
            
            delta1 = sqrt(gamma) * delta;
            
            bb = sqrt(w_2(k, n)) .* [x_tau' ,  - X_QR'];
            
            B = [diag(delta1),b,b1;bb];
                
            for i = 1:Lg*LM
                alpah1 = alpah * conj(B(Lg*LM+1,i)); 
                delta(i) = delta1(i) + alpah1 * B(Lg*LM+1,i);
        
                sigema = delta1(i) / delta(i);
                rho(i) = alpah1 / delta(i);
        
                alpah = alpah * sigema;
                tau0 = B(i,Lg*LM+1);
                tau1 = B(i,Lg*LM+2);
        
                B(i,Lg*LM+1) = sigema * B(i,Lg*LM+1) + rho(i) * B(Lg*LM+1,Lg*LM+1);
                B(i,Lg*LM+2) = sigema * B(i,Lg*LM+2) + rho(i) * B(Lg*LM+1,Lg*LM+2);
        
                B(Lg*LM+1,Lg*LM+1) = B(Lg*LM+1,Lg*LM+1) - tau0 * B(Lg*LM+1,i);
                B(Lg*LM+1,Lg*LM+2) = B(Lg*LM+1,Lg*LM+2) - tau1 * B(Lg*LM+1,i);
            end

            Gamma(Lg*LM) = 0; theta(Lg*LM) = -B(Lg*LM,Lg*LM+1);
            Gamma1(Lg*LM) = 0;theta1(Lg*LM) = -B(Lg*LM,Lg*LM+2);
    
            for i = Lg*LM - 1 :-1:1
                Gamma(i) = Gamma(i+1) + B(Lg*LM+1,i+1) * theta(i+1);
                theta(i) = - B(i,Lg*LM+1) - rho(i) * Gamma(i);

                Gamma1(i) = Gamma1(i+1) + B(Lg*LM+1,i+1) * theta1(i+1);
                theta1(i) = - B(i,Lg*LM+2) - rho(i) * Gamma1(i);        
            end
            G_hat = [theta,theta1];

            
%             D(LM*Lg+1 , :) = sqrt(w_2(k, n)) .* [x_tau' , X_QR'];
% 
%             R_qr = qrgv_1(D, Lg*LM);%QR分解
%             R_qr_MLg =  R_qr(1: Lg*LM, 1: Lg*LM);
%             d_qr_M = R_qr(1:Lg*LM , Lg*LM+1:Lg*LM+LM);
%             for chan = 1:LM
%                 G_hat(:,chan) = backsolution(R_qr_MLg , d_qr_M(:,chan));
%             end
            
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
            
%             D = sqrt(gamma) * R_qr; 
            %********************
        
            %***存储中间变量*****
%             ye=n+(k-1)*LN;
%             parm.G{k,n} = G_hat;
%             parm.R_qr_MLg{k,n} = R_qr_MLg;
%             parm.d_qr_M{k,n} = d_qr_M;
            
            %*******************
        
            %***代价函数********
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


