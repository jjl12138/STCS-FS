function [ lambdaEM, thetaEM, phiEM, lambda_in, p01_upd, lambdat_upd] = CBG_EM( rhat, rvar1, parSig, optEM)
%Formula in Expectation-Maximization Bernoulli-Gaussian Approximate Message
%Passing by Jeremy Vila 

lambda1 = parSig.lambda;%Sparsity
theta1 = parSig.theta;%Mean of the Gaussian distribution
phi1 = parSig.phi;%Variance of the Gaussian distribution
%N = parSig.N;%size of the vector
%T = parSig.T;%size of the vector
[N, T] = size(rhat);
lambda_in1 = parSig.lambda_in;
p01 = parSig.p01;
if strcmp(optEM.struct, 'MC_BitSup')
    lambdat = parSig.lambda_t;
end

lambdaEM = lambda1;
thetaEM = theta1;
phiEM = phi1;
lambda_in = lambda_in1;
p01_upd = p01;
if strcmp(optEM.struct, 'MC_BitSup')
    lambdat_upd = lambdat;
end

lambda = lambda1;
theta = repmat(theta1,N,1);
phi = repmat(phi1,N,1);
rvar = repmat(rvar1,N,1);

% alpha = -abs(rhat - theta).^2 ./ (phi + rvar) + abs(rhat).^2 ./ rvar;
% % bound alpha to avoid overflow 
% alpha = min(max(alpha, -50), 50);
% 
% % for real system
% %beta = lambda./(1-lambda) .* sqrt(rvar./(phi + rvar)) .* exp(alpha);
% % for complex system
% %beta = lambda./(1-lambda) .* rvar./(phi + rvar) .* exp(alpha);
% %pii = 1./(1 + 1./beta);
% 
% % use inbeta relpace 1./beta
% inbeta = (1-lambda)./lambda.*(phi+rvar)./rvar.*exp(-alpha);
% %inbeta = min(max(inbeta, 1e-100), 1e100);
% 
% pii = 1./(1+inbeta);
% %pii(isnan(pii)) = 0.999;
% 
% gamma = (rhat./rvar + theta./phi) ./ (1./rvar + 1./phi);
% nu = 1./(1./rvar + 1./phi);




alpha = - abs(rhat - theta).^2 ./ (phi + rvar) + abs(rhat).^2 ./ rvar;
alpha = min(max(alpha, -100), 100);
%alphaObs = exp(alpha);
beta = (phi+rvar)./rvar.*exp(-alpha);
beta = 1./(1+beta);

switch optEM.struct
    
    case 'MC'
        pi_out = beta;
        [pi_in, s_post, p01_upd] = SPD_MC(pi_out, lambda1(1,1), p01);
        %pi_in should give the estimator
        %p01_upd should be saved

        pii = s_post;
        lambda_upd = sum(pii, 1) / N;
        %can not use s_post(1,1), since lambda_upd is used to calculate p10
        %in Markov Chain process
        lambda_upd = min(max(eps, lambda_upd), 1);
        lambda_upd = ones(N,1) * lambda_upd;
        pii_tmp = pii ./ lambda_upd;

        lambda_in = pi_in;
        
    case 'MC_MMV'
        pi_out = beta;
        [pi_in, s_post, p01_upd] = SPD_MC_MMV(pi_out, mean(lambda1), p01);
        %pi_in should give the estimator
        %p01_upd should be saved

        pii = s_post;
        lambda_upd = sum(pii, 1) / N;
        lambda_upd = min(max(0, lambda_upd), 1);
        lambda_upd = ones(N,1) * lambda_upd;
        pii_tmp = pii ./ lambda_upd;

        lambda_in = pi_in;
      
    case 'MC_ComSup'
        pi_out = beta;
        [pi_in, s_post, p01_upd] = SPD_MC_ComSup1(pi_out, mean(lambda1), p01);
        %pi_in should give the estimator
        %p01_upd should be saved

        pii = s_post;
        lambda_upd = sum(pii, 1) / N;
        lambda_upd = min(max(0, lambda_upd), 1);
        lambda_upd = ones(N,1) * lambda_upd;
        pii_tmp = pii ./ lambda_upd;

        lambda_in = pi_in;
    
    case 'MC_BitSup'
        pi_out = beta;
        [pi_in, s_post, p01_upd, lambdat_upd] = SPD_MC_BitSup(pi_out, mean(lambda1), p01, lambdat);
        %pi_in should give the estimator
        %p01_upd should be saved

        pii = s_post;
        lambda_upd = sum(pii, 1) / N;
        lambda_upd = min(max(0, lambda_upd), 1);
        lambda_upd = ones(N,1) * lambda_upd;
        pii_tmp = pii ./ lambda_upd;

        lambda_in = pi_in;
        
    case 'NoStruct'
        pii = (lambda .* beta) ./ (lambda .* beta + ...
                                (1 - lambda) .* (1 - beta));
        lambda_upd = sum(pii)/N;
        lambda_upd = min(max(0, lambda_upd), 1);
        lambda_upd = ones(N,1) * lambda_upd;

        pii_tmp = pii ./ lambda_upd;
    
    otherwise
        error('Unrecognized Signal Structure');
    
end
                    
gamma = (rhat./rvar + theta./phi) ./ (1./rvar + 1./phi);
nu = 1./(1./rvar + 1./phi);



% Sparsity updata
if optEM.sparsity
    switch optEM.struct     
        case 'MC'
            %lambdaEM = ones(N,1) * lambda_upd;
            lambdaEM = lambda_upd;
        case 'MC_MMV'
            lambdaEM = lambda_upd;
        case 'MC_ComSup'
            lambdaEM = lambda_upd;
        case 'MC_BitSup'
            lambdaEM = lambda_upd;
        case 'NoStruct'
            %lambdaEM = ones(N,1) * lambda_upd ;
            lambdaEM = lambda_upd ;
            lambda_in = lambdaEM;
            
        otherwise
            error('Unrecognized Signal Structure');
            
    end 
    
end

% Mean update
if optEM.mean
    thetaEM = sum(pii_tmp.*gamma)/N;
end

% Variance update
if optEM.var
    switch optEM.struct
        case 'MC_ComSup'
            phiEMtmp = sum(sum( pii_tmp .* ( abs(gamma - theta).^2 + nu ) ))/N/T;     
            phiEM = repmat(phiEMtmp,1,T);
        otherwise
            phiEM = sum( pii_tmp .* ( abs(gamma - theta).^2 + nu ) )/N;
            %phiEM = min(max(phiEM, 0), 10);  
    end
end



end