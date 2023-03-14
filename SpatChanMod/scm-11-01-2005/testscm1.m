clear;
clc;
rand('state',112358);
randn('state',112358);

M=128;
FC=dftmtx(M)/sqrt(M);
MaxiterH=200;

%H=scm(scmparset,linkparset(1),antparset);
scmpar=scmparset;
scmpar.NumBsElements=M;
scmpar.NumMsElements=1;
scmpar.NumPaths=6;%5
scmpar.NumTimeSamples=MaxiterH;
scmpar.Scenario='suburban_macro';
%scmpar.Scenario='urban_macro';
%scmpar.BsUrbanMacroAS='fifteen';
%scmpar.BsUrbanMacroAS='eight';
linkpar = linkparset(1);
%linkpar.MsBsDistance = 50;
antpar = antparset;
%linkpar.MsVelocity=100;
H1=scm(scmpar,linkpar,antpar);
for itert=1:MaxiterH
%    Htem = H1(1,:,1,itert,1);
    Htem = H1(1,:,:,itert,1);
    Htem = sum(Htem,3);
    varH = var(Htem);
    Htem = Htem/sqrt(varH);
    %var(Htem)
    %Htlt(1,:,itert)=Htem;
    Hwtlt(1,:,itert)=Htem*FC;
end;

% for m=1:M
%     Ast(m)=std(Hwtlt(1,m,:))^2;
% end;

% th=mean(Ast(find(Ast>0.1*max(Ast))));
% supp=find(Ast>0.1*th);
% N=length(supp);
% Astsb=zeros(1,M);
% Astsb(supp)=Ast(supp);

hold on;
for iter = 1:MaxiterH
plot(abs(Hwtlt(1,:,iter)))
end
%plot(abs(Hwtlt(1,:,200)))
