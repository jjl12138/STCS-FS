%% Structured Turbo Compressed Sensing algorithm
%   Base on Turbo Compressed Sensing with Partial DFT Sensing Matrix 
%   by Junjie Ma
%   This code is only for partial orthogonal sensing matrix
%%
clc;
clear;

basePath = [fileparts(mfilename('fullpath')) filesep];
addpath([basePath '../pubFun']) %add public function

% N = 128;
% M = 51;
% K = 32;
N = 2000;          % 基站端（发送端）端N个天线（信号长度）
M = 800;           % M个训练时隙
K = 500;           % K个非零值

rand('state',123);  % 均匀分布
randn('state',123); % 高斯分布

Ite_Max = 100;     % 最大迭代次数：100次
NSIM = 10;         % 仿真次数
SNR = 10;          % 信噪比：10dB

NMSE_bp_mea = zeros(Ite_Max, NSIM, length(SNR));   % 每次迭代的NMSE
algTim = zeros(NSIM, length(SNR));                 % 每次迭代耗费的时间

optEM.sparsity = false;
optEM.mean = false;
optEM.var = false;

MC_Supp = true;          % 是否是马尔科夫随机过程
    
for snri = 1:length(SNR)
    
for nsim=1:NSIM  
                                    % 信号的有关参数
    sigPar.theta = 0;               % 高斯分布的均值
    sigPar.phi = N/K;               % 高斯分布的方差
    sigPar.K = K;                   % 导频子载波的个数
    sigPar.N = N;                   % BS端2000个天线
    sigPar.grpNum = 2;              % 非零块的个数
    sigPar.T = 1;
    sigPar.lambda = K/N;            % 稀疏度
    sigPar.lambda_in = K/N;         % ？
    p01 = 0.0040;% for N = 2000, K = 500, grpNum = 2     转移概率
    %p01 = 0.0625;% for N = 128, K = 32, grpNum = 2
    %p01 = 0.03125;% for N = 128, K = 32, grpNum = 1
 
    x = BG_Gen( sigPar, 'CBlock');           % 块稀疏的复伯努利高斯分布
    index = find(x~=0);                      % 找出x中不等于0的值得索引
    A = senMat_Gen(M, N, 'pDFT_RP');         % 生成传感矩阵
    z = A * x;
    if z==0
        nuw = 10^(-SNR(snri)/10);
    else
        nuw = norm(z(:))^2/M*10^(-SNR(snri)/10);
    end
    
    y = z + sqrt(nuw)*(randn(size(z))+1i*randn(size(z)))/sqrt(2);
    
    % record the time strat
    tstart = tic;
    
    vA_pri = 1;
    xA_pri = zeros(N,1);
    
    for it = 1:Ite_Max  
        % LMMSE estimator
        xA_post = xA_pri + vA_pri / ( vA_pri + nuw ) * ( A' * ( y - A * xA_pri) );
        vA_post = vA_pri - M/N * vA_pri^2 / (vA_pri + nuw);
        
        % Update extrinsic messages
        %vA2B_ext = N/M * ( vA_pri + nuw) - vA_pri;
        vA2B_ext = 1/( 1/vA_post - 1/vA_pri );
        vA2B_ext = min(1e11, max(1e-11, vA2B_ext));    % vA2B_ext在[1e-11,1e11]中时才选它，vA2B_ext>1e11时，vA2B_ext=1e11；vA2B_ext<1e-11时，vA2B_ext=1e-11
        xA2B_ext = vA2B_ext/vA_post*xA_post - vA2B_ext/vA_pri*xA_pri;
        
        xB_pri = xA2B_ext;
        vB_pri = vA2B_ext;
        
        % MMSE estimator
        %rhat = xB_pri;
        %rvar = vB_pri;
        rhat = xA2B_ext;
        rvar = vA2B_ext;
        %if  MC_Supp && it>1 % In the case N = 128 
        if  MC_Supp % In the case N = 2000
            %pi_out = Supp_Est( rhat, rvar, sigPar);
            pi_out = CSupp_Est( rhat, rvar, sigPar); 
            %[pi_in, ~, p01_upd] = SPD_MC(pi_out, mean(sigPar.lambda), p01);
            [pi_in, ~, p01_upd] = SPD_MC(pi_out, K/N, p01);
            %[pi_in, ~, p01_upd] = SPD_MC_Cir(pi_out, K/N, p01);
            
            %p01 = p01_upd
            %p01 = 0.1
            %p01 = min(max(p01, 0.001), 0.5);
            sigPar.lambda_in = pi_in; 
        end
        
        [ xhat, xvar, ~ ] = CBG_Est( rhat, rvar, sigPar);
        xB_post = xhat;
        vB_post = mean(xvar);
        %vB_post = xvar;
        
        % Update extrinsic messages
        vB2A_ext = 1/(1/vB_post - 1/vB_pri);
        vB2A_ext = min(1e11, max(1e-11, vB2A_ext));
        xB2A_ext = vB2A_ext/vB_post*xB_post - vB2A_ext/vB_pri*xB_pri;        
        vA_pri = vB2A_ext;
        xA_pri = xB2A_ext;
        
        NMSE_bp_mea(it,nsim,snri) = (norm(x - xhat,'fro')/norm(x,'fro'))^2;
        %fprintf('it = %d\n',it);
    end
    
    % record the time end
    algTim(nsim,snri) = toc(tstart);
    
    fprintf('==================\n');
    fprintf('NSIM:%d, SNR:%ddB\n',nsim,SNR(snri));
    fprintf('Iteration end.\n');
    fprintf('NMSE = %f = %f dB\n',NMSE_bp_mea(Ite_Max,nsim,snri),10*log10(NMSE_bp_mea(Ite_Max,nsim,snri)));
    fprintf('Time = %f\n',algTim(nsim,snri));

end

end

%semilogy(mean(NMSE_bp_mea,2));
plot(10*log10(mean(NMSE_bp_mea,2)));
hold on;
