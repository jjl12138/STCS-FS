function [ pi_out ] = CSupp_Est( rhat, rvar1, parSig)
%CSupp_Est 此处显示有关此函数的摘要
%   此处显示详细说明
%Formula (13) with lambda = 0.5 in Expectation-Maximization Bernoulli-Gaussian Approximate Message
%Passing by Jeremy Vila 
%coded by Chen Lei

N = parSig.N;%size of the vector
%lambda = parSig.lambda;%Sparsity
theta = repmat(parSig.theta, N, 1);%Mean of the Gaussian distribution
phi = repmat(parSig.phi, N, 1);%Variance of the Gaussian distribution
%N = parSig.N;%size of the vector
rvar = repmat(rvar1, N, 1);

tmpVar = rvar + phi;
%pi_out = sqrt(tmpVar./rvar) .* exp(0.5 * ((rhat - theta).^2 ./tmpVar - rhat.^2 ./rvar ));
pi_out = tmpVar./rvar .* exp( abs(rhat - theta).^2 ./tmpVar - abs(rhat).^2 ./rvar );
%0.5* is removed for complex estimation
pi_out = 1./ (1 + pi_out);
end

