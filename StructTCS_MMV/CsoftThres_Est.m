function [ f, fprime ] = CsoftThres_Est( theta, c, ~ )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
%   Asymptotic Analysis of Complex LASSO via Complex Approximate Message
%   Passing Equ.(4)
[N, T] = size(theta);

f = zeros(N,T);
fprime = 0;
thetaAmp = abs(theta);
thetaAng = angle(theta);

index1 = thetaAmp>c;

%f(index1) = (thetaAmp(index1)-c).*cos(thetaAng(index1)) + 1i*(thetaAmp(index1)-c).*sin(thetaAng(index1));
f(index1) = theta(index1) - c./thetaAmp(index1).*theta(index1);
end

