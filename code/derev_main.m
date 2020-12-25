%adaptive-MCLP-dereverberation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% METHODS 
% 1       RLS-based MCLP [1]
% 2       QRRLS-based MCLP [2]
% 3       RLS-based MCLP with time-varying forgetting factor [1,3]
% 4       QRRLS-based MCLP with time-varying forgetting factor [6]
% 5       constrained RLS-based MCLP [1]
% 6       constrained QRRLS-based MCLP [1,2]
% 7       TAQRRLS-based MCLP [1,5]

% REFERENCE
% [1] Jukić A, van Waterschoot T, Doclo S. Adaptive speech dereverberation using constrained sparse 
%     multichannel linear prediction[J]. IEEE Signal Processing Letters, 2016, 24(1): 101-105.
% [2] Sun X, Zhou Y, Shu X. Multi-Channel Linear Prediction Speech Dereverberation Algorithm Based 
%     on QR-RLS Adaptive Filter[C]//Proceedings of the 3rd International Conference on Multimedia 
%     Systems and Signal Processing. 2018: 109-113.
% [3] Zhou Y, Chan S C, Ho K L. A new variable forgetting factor QR-based recursive least M-estimate 
%     algorithm for robust adaptive filtering in impulsive noise environment[C]//2006 14th European 
%     Signal Processing Conference. IEEE, 2006: 1-5.
% [4] Chan S C, Yang X X. Improved approximate QR-LS algorithms for adaptive filtering[J]. IEEE Transactions 
%     on Circuits and Systems II: Express Briefs, 2004, 51(1): 29-39.
% [5] Xinyu Tang, Yang Xu, Rilin Chen and Yi Zhou. A Time-Varying Forgetting Factor-Based QRRLS Algorithm for 
%     Multichannel Speech Dereverberation. 2020 IEEE International Symposium on Signal Processing and 
%     Information Technology (ISSPIT). IEEE, 2020.

% MODIFICATIONS
%      time         author        description
% 1.   2020.12.25   Xinyu Tang    initial version
% 2.   
% 3.   
%% 
clear all;
close all;
%% INPUT
input_dir = 'input\';
output_dir = 'output\';

dirOutput = dir(fullfile(input_dir,'*.wav'));
file_name = {dirOutput.name};
file_name = cellfun(@(c)c(1:end-4) ,file_name,'UniformOutput',false).';
% file_name = {'dev2_1_2_3_RT1.2'};
%% PARAM
wlen = 512;%1024; %fft
gamma =0.94; %forgetting factor
Lg = 4;  %length of filter
T60 = 0.5; %RT60 (s)

gamma_lower = 0.96; %lower bound of forgetting factor
gamma_upper = 0.999; %upper bound of forgetting factor

%% DEREV-METHOD
for i = 1:size(file_name,1)
    %audioread
    [x,fs]=audioread([input_dir,file_name{i},'.wav']);
    mkdir([output_dir file_name{i}]);
    
    for method = 4
        if method == 1
            loop = 1;
            for gamma = [0.99]; %forgetting factor
                [y,e_rls,e_w_rls,J_rls,parm_rls] = derev_MCLP_RLS(x,wlen,gamma,Lg,T60);
                audiowrite([output_dir,'\',file_name{i},'\',file_name{i},'_derevRLS_gamma',num2str(gamma),'_Lg',num2str(Lg),'_T60',num2str(T60),'(10^-4).wav'], y, fs);
            end
        elseif method == 2
            loop = 1;
            for gamma =[0.96 0.999]; %forgetting factor 0.999
                [y,e_qr,e_w_qr,J_qr,parm_qr] = derev_MCLP_QRRLS(x,wlen,gamma,Lg,T60);
                audiowrite([output_dir,'\',file_name{i},'\',file_name{i},'_derevQRRLS_gamma',num2str(gamma),'_Lg',num2str(Lg),'_T60',num2str(T60),'(10^2)l.wav'], y, fs);
            end
        elseif method == 3
            [y,gam]  = derev_MCLP_VFFRLS(x,wlen,gamma_lower,gamma_upper,Lg,T60);
            audiowrite([output_dir,'\',file_name{i},'\',file_name{i},'_derevVFFRLS_gamma',num2str(gamma_lower),num2str(gamma_upper),'_Lg',num2str(Lg),'_T60',num2str(T60),'.wav'], y, fs);
        elseif method == 4
            [y,gam,e_VFF,e_w_VFF,J_VFF,parm_VFF]  = derev_MCLP_VFFQRRLS(x,wlen,gamma_lower,gamma_upper,Lg,T60);
            audiowrite([output_dir,'\',file_name{i},'\',file_name{i},'_derevVFFQRRLS_gamma',num2str(gamma_lower),num2str(gamma_upper),'_Lg',num2str(Lg),'_T60',num2str(T60),'(10^2).wav'], y, fs);            
        elseif method == 5
            for gamma = [0.96 0.999];
                [y_rls_u,y_rls_z] = derev_MCLP_RLS_ADMM(x,wlen,gamma,Lg,T60);
                audiowrite([output_dir,'\',file_name{i},'\',file_name{i},'_derevADMMRLS_u_gamma',num2str(gamma),'_Lg',num2str(Lg),'_T60',num2str(T60),'(10^-4).wav'], y_rls_u, fs);
                audiowrite([output_dir,'\',file_name{i},'\',file_name{i},'_derevADMMRLS_z_gamma',num2str(gamma),'_Lg',num2str(Lg),'_T60',num2str(T60),'(10^-4).wav'], y_rls_z, fs);
            end
        elseif method == 6
            for gamma = [0.94 0.975 0.99];
                [y_qr_u,y_qr_z] = derev_MCLP_QRRLS_ADMM(x,wlen,gamma,Lg,T60);
                audiowrite([output_dir,'\',file_name{i},'\',file_name{i},'_derevADMMQRRLS_u_gamma',num2str(gamma),'_Lg',num2str(Lg),'_T60',num2str(T60),'.wav'], y_qr_u, fs);
                audiowrite([output_dir,'\',file_name{i},'\',file_name{i},'_derevADMMQRRLS_z_gamma',num2str(gamma),'_Lg',num2str(Lg),'_T60',num2str(T60),'.wav'], y_qr_z, fs);
            end
        elseif method == 7
            gamma = 0.96;
            [y,e_qr,e_w_qr,J_qr,parm_qr] = derev_MCLP_TAQRRLS(x,wlen,gamma,Lg,T60);
            audiowrite([output_dir,'\',file_name{i},'\',file_name{i},'_derevTAQRRLS_gamma',num2str(gamma),'_Lg',num2str(Lg),'_T60',num2str(T60),'(10^2).wav'], y , fs);
        end  
    end
 end
