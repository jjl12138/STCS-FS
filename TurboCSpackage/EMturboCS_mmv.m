%%  STCS combinating with Expectation Maximization
%   also combine with different delay domain denoiser
%%
clc;
clear;

basePath = [fileparts(mfilename('fullpath')) filesep];
addpath([basePath '../pubFun']); %add public function

rand('state',123);
randn('state',123);

N = 256;
K = 64;
T = 32;

%SNR = 14;
% kSIM = 0.98;
% Mpt = 128;
SNR = 30;
%kSIM = 0.05:0.05:1;
% SNR = 10:2:30;
kSIM = 0.4;
Mpt = ceil(kSIM*N);

% SNR = 60;
% %kSIM = 0.02:0.02:1;
% % Mpt = ceil(kSIM*N);
% Mpt = 1:1:256;

% load Hp_Ang_Fre1_100.mat;
% load Hp_Ang_Del1_100.mat;

NSIM = 10;
%NMSE = nan(NSIM,1);
NMSEreChan = zeros(NSIM, length(Mpt), length(SNR));

for mi = 1:length(Mpt)

for snri = 1:length(SNR)
    
for nsim = 1:NSIM
    
    optEM.sparsity = true;
    optEM.mean = false;
    optEM.var = true;
    %optEM.struct = 'NoStruct';
    optEM.struct = 'MC_ComSup';
    %optEM.struct = 'MC_BitSup';
    
    %x = cell(1,T);
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
    Xt = [xt{:}];
    X = [x{:}];

    %X = HAngFre(:,:,nsim);
    %X = X/sqrt(var(X(:)));
    x = cell(1,T);
    xt = cell(1,T);
    for t = 1:T
       x{t} = X(:,t);
       xt{t} = Xt(:,t);
    end

    M = Mpt(mi);
    for t = 1:T
    %         x{t} = BG_Gen( sigPar, 'CBlock');
    %         index{t} = find(x{t}~=0);

        %this code is only for partial orthogonal sensing matrix
        A{t} = senMat_Gen(M, N, 'pDFT_RP');   
        if strcmp(optEM.struct, 'MC_BitSup')
            z{t} = A{t} * xt{t};
        else
            z{t} = A{t} * x{t};
        end
    end

    ztmp = [z{:}];
    varw = norm(ztmp(:))^2/M/T;

    for t = 1:T
        nuw{t} = varw*10^(-SNR(snri)/10);%ÐèÒªÐÞ¸Änuw = norm(u(:))^2/M*10^(-isnr/10);
        %y = z + sqrt(nuw)*randn(size(z));
        y{t} = z{t} + sqrt(nuw{t})*(randn(size(z{t}))+1i*randn(size(z{t})))/sqrt(2);
    end

    %end

    sigEst.theta = 0*ones(1,T);
    %sigEst.lambda = (M/2)/N*ones(N,T);
    %sigEst.lambda_in = (M/2)/N*ones(N,T);
    sigEst.lambda = 0.3*ones(N,T);
    sigEst.lambda_in = 0.3*ones(N,T);
    sigEst.p01 = 0.1*ones(1,T);
    if strcmp(optEM.struct, 'MC_BitSup')
        sigEst.lambda_t = 0.1*ones(1,T);
    end
    sigEst.N = N;
    for t = 1:T
        sigEst.psi(1,t) = norm([y{t}])^2 / ((100+1)*M);
        sigEst.phi(1,t) = norm([y{t}])^2 / ( norm([A{t}],'fro')^2*((M/2)/N) );
        %if strcmp(optEM.struct, 'MC_BitSup')
        %    sigEst.phi(1,t) = norm([y{t}])^2 / ( norm([A{t}],'fro')^2*((M/2)/N) );
        %end
    end   


    vA_pri = 1*ones(1,T);
    for t = 1:T
        vA_pri(1,t) = norm(y{t},'fro')^2/M; 
        %if strcmp(optEM.struct, 'MC_BitSup')
        %    vA_pri(1,t) = norm(y{t},'fro')^2/M; 
        %end
    end
    xA_pri = zeros(N,T);

    %if snri==1 && mi==1 && nsim==42
    %sigEst = sigPar;

    smooth_iters = 5;%5
    opt.nit = 30;%30
    opt.timSupp = false;
    opt.step = 1;
    %opt.nuw = sigEst.psi;
    opt.nuw = nuw;
    
%     if SNR(snri)>12
%     opt.nuw{1} = varw*10^(-12/10);
%     opt.nuw{2} = varw*10^(-12/10);
%     opt.nuw{9} = varw*10^(-14/10);
%     opt.nuw{11} = varw*10^(-14/10);
%     opt.nuw{13} = varw*10^(-14/10);
%     %opt.nuw{52} = varw*10^(-20/10);
%     end
%     if SNR(snri)>16
%     opt.nuw{1} = varw*10^(-10/10);
%     opt.nuw{2} = varw*10^(-10/10);
%     opt.nuw{9} = varw*10^(-10/10);
%     opt.nuw{11} = varw*10^(-10/10);
%     opt.nuw{13} = varw*10^(-10/10);
%     %opt.nuw{52} = varw*10^(-20/10);
%     end
    
    opt.xA_pri0 = [];
    opt.vA_pri0 = [];
    opt.vA_pri0 = vA_pri;
    
    opt.saveHist = false;
    if opt.saveHist == true
        warning('WARNING: Data is SAVING!!');
    end
    opt.x = x;
    opt.tol = 1e-20;

    %[opt.U,opt.D] = eig(A*A');
    %opt.d = diag(opt.D);
    opt.U = []; opt.D = []; opt.d = [];

    estHist2.xA_post = [];
    estHist2.vA_post = [];
    estHist2.xB_pri = [];
    estHist2.vB_pri = [];
    estHist2.xB_post = [];
    estHist2.vB_post = [];
    estHist2.xA_pri = [];
    estHist2.vA_pri = [];
    estHist2.nuw = [];
    estHist2.errX = [];
    estHist2.difX = [];
    estHist2.difX2 = [];

    estHist2.lambda = [];
    estHist2.theta = [];
    estHist2.phi = [];
    estHist2.lambda_in = [];
    estHist2.p01 = [];

    stop = false;
    i = 0;
    while ~stop

        % Increment time and check if EM iterations are done
        i = i + 1;
        if i >= smooth_iters
            stop = true;
        end
        
        [ estFin, optFin, estHist] = tcsEst_mmv( y, A, opt, sigEst );

        if strcmp(optEM.struct, 'MC_BitSup')
            [ sigEst.lambda, sigEst.theta, sigEst.phi, sigEst.lambda_in, sigEst.p01, sigEst.lambda_t ] ...
                = CBG_EM( estFin.xB_pri, estFin.vB_pri, sigEst, optEM);
            %opt.xA_pri0 = estFin.xA_pri;
            %opt.vA_pri0 = estFin.vA_pri;
        else
            [ sigEst.lambda, sigEst.theta, sigEst.phi, sigEst.lambda_in, sigEst.p01 ] ...
                = CBG_EM( estFin.xB_pri, estFin.vB_pri, sigEst, optEM);
            
        end
        %opt.xA_pri0 = estFin.xA_pri;
        %opt.vA_pri0 = estFin.vA_pri;
        %opt.nuw = estFin.nuw;
        xhat = estFin.xB_post;

        if opt.saveHist 
            estHist2.xA_post = [estHist2.xA_post, estHist.xA_post];
            estHist2.vA_post = [estHist2.vA_post, estHist.vA_post];
            estHist2.xB_pri = [estHist2.xB_pri, estHist.xB_pri];
            estHist2.vB_pri = [estHist2.vB_pri, estHist.vB_pri];
            estHist2.xB_post = [estHist2.xB_post, estHist.xB_post];
            estHist2.vB_post = [estHist2.vB_post, estHist.vB_post];
            estHist2.xA_pri = [estHist2.xA_pri, estHist.xA_pri];
            estHist2.vA_pri = [estHist2.vA_pri, estHist.vA_pri];
            estHist2.nuw = [estHist2.nuw, estHist.nuw];
            estHist2.errX = [estHist2.errX, estHist.errX];
            estHist2.difX = [estHist2.difX, estHist.difX];
            estHist2.difX2 = [estHist2.difX2, estHist.difX2];

            estHist2.lambda = [estHist2.lambda, sigEst.lambda];
            estHist2.theta = [estHist2.theta, sigEst.theta];
            estHist2.phi = [estHist2.phi, sigEst.phi];
            estHist2.lambda_in = [estHist2.lambda_in, sigEst.lambda_in];
            estHist2.p01 = [estHist2.p01, sigEst.p01];
        end
    end

    %NMSE(nsim,1) = (norm(x - xhat,'fro')/norm(x,'fro'))^2;
    %NMSE = (norm(x - xhat,'fro')/norm(x,'fro'))^2;
    if strcmp(optEM.struct, 'MC_BitSup')
        NMSE = (norm([xt{:}] - xhat,'fro')/norm([xt{:}],'fro'))^2;
    else
        NMSE = (norm([x{:}] - xhat,'fro')/norm([x{:}],'fro'))^2;
    end
    
    if isnan(NMSE)
        NMSE = 1;
    end
    NMSEreChan(nsim, mi, snri) = NMSE;

    fprintf('SNR = %d dB, M = %d, nsim = %d ',SNR(snri),M,nsim);
    fprintf('NMSE = %.4f dB \n',10*log10(NMSE));
    fprintf('Iteration end.\n');
    

end
% fprintf('SNR = %d dB, M = %d, nsim = %d \n',SNR(snri),M,nsim);
% end

end

end

%save NMSEphaNorea_mmWave NMSEreChan;