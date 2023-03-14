function [ x ] = BG_Gen( parSig, prior)
%genBG 
%The detail goes here.
%coded by Lei Chen

K = parSig.K;  % 稀疏度
%lambda = parSig.lambda;%Sparsity
theta = parSig.theta(1,1);%Mean of the Gaussian distribution
phi = parSig.phi(1,1);%Variance of the Gaussian distribution
N = parSig.N;%size of the vector
T = parSig.T;

switch prior
    % Standard Bernoulli Gaussian distribution
    case 'BG'
        %index = random('Binomial',1,K/N,N,1);
        loc = randperm(N);
        loc = loc(1,1:K);
        index = zeros(N,1);
        index(loc) = 1;
        x = sqrt(phi) * randn(N,1) + theta;
        x = x.*index;
        
    % Complex Bernoulli Gaussian distribution
    case 'CBG'
        %index = random('Binomial',1,K/N,N,1);
        x = zeros(N,T);
        for i = 1:T
            loc = randperm(N);
            loc = loc(1,1:K(1,i));
            index = zeros(N,1);
            index(loc) = 1;
            x(:,i) = sqrt(phi(1,i)) * (randn(N,1) + 1i*randn(N,1))/sqrt(2) + theta(1,i);
            x(:,i) = x(:,i).*index;
        end
        
    % Bernoulli Gaussian distribution with Block structure
    case 'Block'
        %K = fix(lambda*N);      % sparsity (nonzero elements in the signal)
        q = parSig.grpNum;       % block number (number of nonzero blocks)
        
        % generate sparse signal
        x=zeros(N,1);
        r=abs(randn(q,1));% generate non-negative random number 
        r=r+1; r=round(r*K/sum(r));% divide K into q blocks 
        r(q)=K-sum(r(1:q-1));% revise the last block and make sure the number of nonzeros is K 
        g=round(r*N/K);% locate the blocks
        g(q)=N-sum(g(1:q-1));% make sure the location is reasonable
        g_cum=cumsum(g);% give the right location

        for ii=1:q,

            % intra-block correlation
            %beta = 1 - 0.05*rand;  
            beta = 0;

            % generate i-th group with intra-block correlation
            seg = [];
            seg(1) = randn;
            for kk = 2 : r(ii)
                seg(kk,1) = beta*seg(kk-1,1) + sqrt(1-beta^2) * randn;% a Markov Random Process
            end
            loc=randperm(g(ii)-r(ii)-2);
            x_tmp=zeros(g(ii), 1);
            x_tmp(loc(1)+2:loc(1)+1+r(ii))= seg; % generate coherent blocks

            x(g_cum(ii)-g(ii)+1:g_cum(ii), 1)=x_tmp;
        end
        x = sqrt(phi) * x + theta;
        
    % Complex  Bernoulli Gaussian distribution with Block structure
    case 'CBlock'
        %K = fix(lambda*N);      % sparsity (nonzero elements in the signal)
        q = parSig.grpNum;       % group number (number of nonzero blocks)
        
        % generate sparse signal
        x=zeros(N,1) ;
        r=abs(randn(q,1));% generate non-negative random number 
        r=r+1; 
        r=round(r*K/sum(r));% divide K into q blocks 
        r(q)=K-sum(r(1:q-1));% revise the last block and make sure the number of nonzeros is K 
        g=round(r*N/K);% locate the blocks
        g(q)=N-sum(g(1:q-1));% make sure the location is reasonable
        g_cum=cumsum(g);% give the right location

        for ii=1:q,

            % intra-block correlation
            beta = 1 - 0.05*rand;  
            %beta = 0;

            % generate i-th group with intra-block correlation
            seg = [];
            seg(1) = sqrt(phi) * (randn + 1i*randn)/sqrt(2);
            for kk = 2 : r(ii)
                seg(kk,1) = beta*seg(kk-1,1) + sqrt(1-beta^2) * sqrt(phi) * (randn + 1i*randn)/sqrt(2);% a Markov Random Process
            end
            loc=randperm(g(ii)-r(ii)-2);
            x_tmp=zeros(g(ii), 1);
            x_tmp(loc(1)+2:loc(1)+1+r(ii))= seg; % generate coherent blocks

            x(g_cum(ii)-g(ii)+1:g_cum(ii), 1)=x_tmp;
        end
        x = x + theta;
        
    %Bernoulli Gaussian distribution with group structure
    case 'Grp'
        grpNum = parSig.q;
        grpSize = K/grpNum;
        perm = randperm(N/grpSize);
        index_sel = perm(1:grpNum);
        tmp = zeros(1,N/grpSize);
        tmp(1,index_sel) = 1;
        tmp = repmat(tmp, grpSize, 1);
        tmp = reshape(tmp, N, 1);
        x = sqrt(phi)*randn(N,1) + theta;
        x = x.*tmp;
      
    otherwise
        error('Unrecognized marginal signal prior');
end

