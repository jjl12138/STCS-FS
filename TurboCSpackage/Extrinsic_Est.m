function [ xA_pri ] = Extrinsic_Est( xB_post, xB_pri, vB_pri, sigEst )
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明

N = size(xB_post,1);
T = size(xB_post,2);
F = dftmtx(T)/sqrt(T);
LTimDom = 16;

delta = 1e-2;%1e-2
inIter = 5;%5
div = NaN*ones(inIter, T);
cOpt1 = NaN*ones(1,T);
xA_pri = NaN*ones(N,T);

for inIt = 1:inIter
Ntil = randn(N,T) + 1i*randn(N,T)/sqrt(2);
[ xhat1, xvar, ~ ] = CBG_Est( xB_pri + delta*Ntil, repmat(vB_pri, N, 1), sigEst);
xvar1 = mean(xvar, 1);
xBt_post =  xhat1*(F').';

% Soft Thresholding Denoiser
[ xBt_post, ~ ] = CsoftThres_Est( xBt_post, sqrt(mean(xvar1,2)), []);
% Cut-off Denoiser
%xBt_post(:,LTimDom+1:1:end) = 0;

xhat1 = xBt_post*F.';

diff = (xhat1 - xB_post)/delta;

div(inIt, :) = sum(diff.*Ntil,1);
end
divMean = mean(div,1);
cOpt2 = divMean/N;

for t = 1:T
    cOpt1(1,t) = real( xB_pri(:,t)'*(xB_post(:,t) - cOpt2(1,t)*xB_pri(:,t)) )/( (xB_post(:,t)-cOpt2(1,t)*xB_pri(:,t))'*(xB_post(:,t)-cOpt2(1,t)*xB_pri(:,t)) );
end

for t = 1:T
    xA_pri(:,t) = cOpt1(1,t)*(xB_post(:,t) - cOpt2(1,t)*xB_pri(:,t));
end

% This is B2A_Cal.m, and extrinsic calculation is not correct.
%     vB2A_ext = 1/(1/vB_post - 1/vB_pri);
%     vB2A_ext = min(1e11, max(1e-11, vB2A_ext));
%     xB2A_ext = vB2A_ext/vB_post*xB_post - vB2A_ext/vB_pri*xB_pri;
% 
%     vA_pri = vB2A_ext;
%     xA_pri = xB2A_ext;
end