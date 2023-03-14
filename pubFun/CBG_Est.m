function [ xhat, xvar, cOpt2 ] = CBG_Est( rhat, rvar, parSig)
%CBG_Est 此处显示有关此函数的摘要
%   此处显示详细说明

N = parSig.N;%size of the vector
%lambda = repmat(parSig.lambda, N, 1);%Sparsity
lambda = parSig.lambda_in;
theta = repmat(parSig.theta, N, 1);%Mean of the Gaussian distribution
phi = repmat(parSig.phi, N, 1);%Variance of the Gaussian distribution

alpha = abs(rhat - theta).^2 ./ (phi + rvar) - abs(rhat).^2 ./ rvar;
% bound alpha to avoid overflow 
%alpha = min(max(alpha, -100), 100);

% for real system
%beta = lambda./(1-lambda) .* sqrt(rvar./(phi + rvar)) .* exp(0.5 * alpha);

% for complex system
% beta = lambda./(1-lambda) .* rvar./(phi + rvar) .* exp(-alpha);%0.5* is removed for complex estiamtion
% pii = 1./(1 + 1./beta);

% use inbeta relpace 1./beta
inbeta = (1-lambda)./lambda.*(phi+rvar)./rvar.*exp(alpha);
%inbeta = min(max(inbeta, 1e-100), 1e100);

pii = 1./(1+inbeta);
%pii(isnan(pii)) = 0.999;

gamma = (rhat./rvar + theta./phi) ./ (1./rvar + 1./phi);
nu = 1./(1./rvar + 1./phi);

xhat = pii.*gamma;
xvar = pii.*(nu + abs(gamma).^2) - pii.^2 .* abs(gamma).^2;

cOpt2 = 0;
%cOpt2 = mean(pii.*phi./(rvar+phi));

% alpha = -abs(rhat - theta).^2 ./ (phi + rvar) + abs(rhat).^2 ./ rvar;
% beta = (phi+rvar)./rvar.*exp(-alpha);
% beta = 1./(1+beta);
% 
% 
% (PI_IN .* PI_OUT) ./ (PI_IN .* PI_OUT + ...
%                         (1 - PI_IN) .* (1 - PI_OUT));

end

