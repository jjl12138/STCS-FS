clc;
clear;

basePath = [fileparts(mfilename('fullpath')) filesep];
addpath([basePath '../pubFun']); %add public function

rand('state',123);
randn('state',123);

N = 2000;
K = 500;

SNR = 10;
Mpt = 800;

NSIM = 5;
%NMSE = nan(NSIM,1);
NMSEreChan = zeros(NSIM, length(Mpt), length(SNR));


for mi = 1:length(Mpt)

for snri = 1:length(SNR)
    
for nsim = 1:NSIM

optEM.sparsity = true;
optEM.mean = false;
optEM.var = true;
%optEM.struct = 'MC_MMV';
optEM.struct = 'MC';
%optEM.struct = 'NoStruct';

sigPar.theta = 0;
sigPar.phi = N/K;
sigPar.K = K;
sigPar.N = N;
sigPar.grpNum = 2;
sigPar.T = 1;
sigPar.lambda = K/N;
sigPar.lambda_in = K/N;

x = BG_Gen( sigPar, 'CBlock');
% index = find(x~=0);

%this code is only for partial orthogonal sensing matrix
M = Mpt(mi);
A = senMat_Gen(M, N, 'pDFT_RP');
%A = senMat_Gen(M, N, 'CGaussN');
z = A * x;
nuw = norm(z(:))^2/M*10^(-SNR(snri)/10);%ÐèÒªÐÞ¸Änuw = norm(u(:))^2/M*10^(-isnr/10);
%y = z + sqrt(nuw)*randn(size(z));
w = sqrt(nuw)*(randn(size(z))+1i*randn(size(z)))/sqrt(2);
y = z + w;


%if snri==1 && mi==1 && nsim==42
%sigEst = sigPar;

sigEst.theta = 0;
sigEst.lambda = 0.3*ones(N,1);%0.3
sigEst.lambda_in = 0.3*ones(N,1);
sigEst.p01 = 0.2;%0.2

% sigEst.theta = 0*ones(1,T);
% sigEst.lambda = (M/2)/N*ones(N,1);
% sigEst.lambda_in = (M/2)/N*ones(N,1);
% sigEst.p01 = 0.1;
%sigEst.lambda = lambda1;
if SNR(snri)<20
    sigEst.psi = norm(y)^2 / ((100+1)*M);
else if SNR(snri)>=20 && SNR(snri)<40
    sigEst.psi = norm(y)^2 / ((1000+1)*M);
    else  
        sigEst.psi = norm(y)^2 / ((10000+1)*M);
    end
end

sigEst.phi = (norm(y)^2 - M*sigEst.psi) / (M * sigEst.lambda(1,1) );
%sigEst.phi = norm([y])^2 / ( norm(A,'fro')^2*((M/2)/N) );
%sigEst.K = K;
sigEst.N = N;

% sigEst.phi = var(x(index));
% sigEst.lambda_in = HmmWave128Spa(nsim,1);
% sigEst.lambda = HmmWave128Spa(nsim,1);

% sigEst.phi = var(x(index));
% sigEst.lambda_in = K/N;
% sigEst.lambda = K/N;
% sigEst.lambda = K/N;



smooth_iters = 5;
opt.nit = 30;
opt.step = 1;
%opt.nuw = sigEst.psi;
opt.nuw = nuw;
% if SNR(snri)>15
%     opt.nuw = norm(z(:))^2/M*10^(-15/10);
% end
opt.xA_pri0 = [];
opt.vA_pri0 = [];
opt.vA_pri0 = norm(y,'fro')^2/M;
opt.saveHist = false;
opt.x = x;
opt.tol = 1e-20;

[opt.U,opt.D] = eig(A*A');
opt.d = diag(opt.D);

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
    
    [ estFin, optFin, estHist ] = tcsEst( y, A, opt, sigEst );
    
    [ sigEst.lambda, sigEst.theta, sigEst.phi, sigEst.lambda_in, sigEst.p01 ] ...
        = CBG_EM( estFin.xB_pri, estFin.vB_pri, sigEst, optEM);
    %opt.xA_pri0 = estFin.xA_pri;
    %opt.vA_pri0 = estFin.vA_pri;
    opt.nuw = estFin.nuw;
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
NMSE = (norm(x - xhat,'fro')/norm(x,'fro'))^2;
NMSEreChan(nsim, mi, snri) = NMSE;

fprintf('SNR = %d dB, M = %d, nsim = %d ',SNR(snri),M,nsim);
fprintf('NMSE = %.4f dB \n',10*log10(NMSE));
fprintf('Iteration end.\n');

end
% fprintf('SNR = %d dB, M = %d, nsim = %d \n',SNR(snri),M,nsim);
% end

end

end







