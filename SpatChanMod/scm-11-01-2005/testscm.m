clear;
clc;

rand('state',1);
randn('state',1);

M=128;
FC=dftmtx(M)/sqrt(M);
MaxiterH=200;

H=scm(scmparset,linkparset(1),antparset);
scmpar=scmparset;
scmpar.NumBsElements=M;
scmpar.NumMsElements=1;
scmpar.NumPaths=6;
scmpar.NumTimeSamples=MaxiterH;
%scmpar.Scenario='urban_macro';
%scmpar.BsUrbanMacroAS='eight';
linkpar=linkparset;
%linkpar.MsVelocity=1;
H1=scm(scmpar,linkpar,antparset);
for itert=1:MaxiterH
    Htem=H1(1,:,1,itert,1);
    %Htlt(1,:,itert)=Htem;
    Hwtlt(1,:,itert)=Htem*FC;
end;

for m=1:M
    Ast(m)=std(Hwtlt(1,m,:))^2;
end;

th=mean(Ast(find(Ast>0.1*max(Ast))));
supp=find(Ast>0.1*th);
N=length(supp);
Astsb=zeros(1,M);
Astsb(supp)=Ast(supp);

