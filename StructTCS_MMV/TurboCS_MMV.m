%% Structured Turbo Compressed Sensing algorithm
%   Developed for Frequency Domain Support Model 
%   Can be combined with delay domain denoiser
%   Based on Turbo Compressed Sensing with Partial DFT Sensing Matrix 
%   by Junjie Ma
%%
clc;
clear;

basePath = [fileparts(mfilename('fullpath')) filesep];
addpath([basePath '../pubFun']) %add public function

N = 256;
M = 103;%205
%M = 128;
K = 64;

rand('state',2);
randn('state',2);

Ite_Max = 50;
NSIM = 10;%1000 
SNR = 30;
T = 32;
F = dftmtx(T)/sqrt(T);

NMSE_bp_mea = zeros(Ite_Max, NSIM, length(SNR));

optEM.sparsity = false;
optEM.mean = false;
optEM.var = false;
optEM.struct = 'MC_ComSup';
%optEM.struct = 'NoStruct';

%MC_Supp = false;
%EM_Supp = false;

TimSupp = false;
    
for snri = 1:length(SNR)
    
for nsim=1:NSIM  
    
    sigPar.theta = zeros(1,T);
    sigPar.phi = N/K*ones(1,T);
    sigPar.K = K;
    sigPar.N = N;
    sigPar.grpNum = 4;
    sigPar.T = T;
    sigPar.lambda = K/N*ones(N,T);
    sigPar.lambda_in = K/N*ones(N,T);
    sigPar.p01 = 0.0625*ones(1,T);
    sigPar.lambda_t = 0*ones(1,T); 
   
    x = cell(1,T);
    index = cell(1,T);
    A = cell(1,T);
    z = cell(1,T);
    nuw = cell(1,T);
    y = cell(1,T);
     
    LTimDom = 16;
    KTimDom = 4;
    KAngDom = K;
    grpNumAngDom = 4;
    [x, xt] = MultiChanGen1( N, T, LTimDom, KTimDom, KAngDom, grpNumAngDom );
    Xt = [xt{:}];%Time domain data matrix
    X = [x{:}];%Frequency domain data matrix
    %[x1,y1] = meshgrid(1:1:T, 1:1:N);
    %mesh(x1,y1, abs(X));

    % MMV signal model
%     xtmp = BG_Gen( sigPar, 'CBlock');
%     index1 = find(xtmp~=0);
%     xtmp = zeros(N,1);
%     for t = 1:T
%         xtmp(index1) = 1;
%         x{t} = sqrt(sigPar.phi(1,t))*(randn(N,1)+1i*randn(N,1))/sqrt(2);
%         x{t} = x{t}.*xtmp;
%     end
%    x{1} = xtmp;
    
    for t = 1:T
        %x{t} = BG_Gen( sigPar, 'CBlock');
        index{t} = find(x{t}~=0);

        %this code is only for partial orthogonal sensing matrix
        %partial DFT can not work well
        %A{t} = senMat_Gen(M, N, 'pInvDFT');
        %A{t} = senMat_Gen(M, N, 'pDFT');
        %A{t} = senMat_Gen(M, N, 'pdouDFT');
        A{t} = senMat_Gen(M, N, 'pDFT_RP');

        z{t} = A{t} * x{t};
        
%         xtmp = x{t};
%         indextmp = index{t};
%         var(xtmp(indextmp))
    end
    
    ztmp = [z{:}];
    varw = norm(ztmp(:))^2/M/T;
    
    for t = 1:T
        nuw{t} = varw*10^(-SNR(snri)/10);%ÐèÒªÐÞ¸Änuw = norm(u(:))^2/M*10^(-isnr/10);
        %y = z + sqrt(nuw)*randn(size(z));
        y{t} = z{t} + sqrt(nuw{t})*(randn(size(z{t}))+1i*randn(size(z{t})))/sqrt(2);
    end
    

    sigEst = sigPar;
%     for t = 1:T
%         sigEst.lambda_in(:,t) = 0.00001*ones(N,1);
%         sigEst.lambda_in(index{t},t) = 1-0.00001;
%     end

    vA_pri = 1*ones(1,T);
    xA_pri = zeros(N,T);
    
    xB_pri = zeros(N,T);
    vB_pri = zeros(1,T);
    
    for it = 1:Ite_Max  
        for t = 1:T
            [ xB_pri(:,t), vB_pri(:,t) ] = A2B_Cal( y{t}, A{t}, nuw{t}, xA_pri(:,t), vA_pri(:,t) );
        end
%         rhat = xB_pri;
%         rvar = vB_pri;

        if it>0
            %[ sigEst.lambda, sigEst.theta, sigEst.phi, sigEst.lambda_in, sigEst.p01 ] ...
            [ ~, ~, ~, sigEst.lambda_in, ~] ...
            = CBG_EM( xB_pri, vB_pri, sigEst, optEM);
        end
        
        [ xhat, xvar, ~ ] = CBG_Est( xB_pri, repmat(vB_pri, N, 1), sigEst);
        xB_post = xhat;
        vB_post = mean(xvar, 1);
        
        if TimSupp
            
            xBt_post = xB_post*(F').';
   
            % Extrinsic Estimator should be revised as well!!!
            % Soft Thresholding Denoiser
            [ xBt_post, ~ ] = CsoftThres_Est( xBt_post, sqrt(mean(vB_post,2)), []);
            % Cut-off Denoiser
            xBt_post(:,LTimDom+1:1:end) = 0;           
            
            xB_post = xBt_post*F.';
            
            [ xA_pri ] = Extrinsic_Est( xB_post, xB_pri, vB_pri, sigEst );
            
            for t = 1:T
                vA_pri(1,t) = ( norm(y{t}-A{t}*xA_pri(:,t))^2 - M*nuw{t})/M;
                vA_pri(1,t) = min(1e11, max(1e-11, vA_pri(1,t)));
            end
                  
        else
            
            for t = 1:T
                [ xA_pri(:,t), vA_pri(1,t) ] = B2A_Cal( xB_post(:,t), vB_post(1,t), xB_pri(:,t), vB_pri(1,t) );
            end
            
        end


        
        % should be changed to xA_pri !!!!
        if TimSupp
            NMSE_bp_mea(it,nsim,snri) = (norm([x{:}] - xA_pri,'fro')/norm([x{:}],'fro'))^2;
        else
            NMSE_bp_mea(it,nsim,snri) = (norm([x{:}] - xhat,'fro')/norm([x{:}],'fro'))^2;
        end
        %fprintf('it = %d\n',it);

    end
    fprintf('==================\n');
    fprintf('NSIM:%d, SNR:%ddB\n',nsim,SNR(snri));
    fprintf('Iteration end.\n');
    fprintf('NMSE = %f, %f dB\n',NMSE_bp_mea(Ite_Max,nsim,snri), 10*log10(NMSE_bp_mea(Ite_Max,nsim,snri)));

end

end

%save TCSComCut10dB NMSE_bp_mea