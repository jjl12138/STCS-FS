function [ PI_IN, S_POST, P01_UPD ] = SPD_MC_ComSup1( PI_OUT, LAMBDA, p01 )
%   SPD_MC_ComSup  
%   forward and backward algorithm
%   binary symbols
%   -1 1
%   0  1
%   This function solves multiple measurement vector with Common Markov
%   Chain support

N = size(PI_OUT,1);
T = size(PI_OUT,2);

P01 = repmat(p01,N,1);
LAMBDA = repmat(LAMBDA,N,1);

P10 = P01 .* (LAMBDA ./ (1 - LAMBDA));
LAMBDA_FWD = NaN*ones(N,1);
LAMBDA_BWD = NaN*ones(N,1);
LAMBDA_FWD(1,:) = LAMBDA(1,1);
LAMBDA_FWD(1,:) = 0.01;
LAMBDA_BWD(N,:) = 0.5;
% LAMBDA_FWD(1,1) = eps;
% LAMBDA_BWD(N,1) = eps;

epsilon = 1e-3;
PI_IN = NaN*ones(N,T);

sumLogPI_OUT = sum(log(PI_OUT),2);
sumLogPI_OUT = max(min(sumLogPI_OUT, 50), -50);
sumLog1MinusPI_OUT = sum(log(1-PI_OUT),2);
sumLog1MinusPI_OUT = max(min(sumLog1MinusPI_OUT, 50), -50);

for n = 2 : 1 : N
    LAMBDA_FWD(n,1) = ( P10(n,1) + (1 - P01(n,1)) * ( LAMBDA_FWD(n-1,1)/(1-LAMBDA_FWD(n-1,1)) ) ...
        * exp(sumLogPI_OUT(n-1,1) - sumLog1MinusPI_OUT(n-1,1)) )/...
        ( ( LAMBDA_FWD(n-1,1)/(1-LAMBDA_FWD(n-1,1)) ) * exp(sumLogPI_OUT(n-1,1) - sumLog1MinusPI_OUT(n-1,1)) + 1 );
end
for n = N-1 : -1 : 1
    LAMBDA_BWD(n,1) = ( P01(n,1) + (1 - P01(n,1)) * ( LAMBDA_BWD(n+1,1)/(1-LAMBDA_BWD(n+1,1)) ) ...
        * exp(sumLogPI_OUT(n+1,1) - sumLog1MinusPI_OUT(n+1,1)) )/...
        ( (1-P01(n,1)+P10(n,1)) * ( LAMBDA_BWD(n+1,1)/(1-LAMBDA_BWD(n+1,1)) ) ...
        * exp(sumLogPI_OUT(n+1,1) - sumLog1MinusPI_OUT(n+1,1)) + (1-P10(n,1)+P01(n,1))  );
end

for n = 1:N
    for t = 1:T
        %tmp = sumLog1MinusPI_OUT(n,1) - log(1-PI_OUT(n,t)) - ( sumLogPI_OUT(n,1) - log(PI_OUT(n,t)) );
        tmp = sumLog1MinusPI_OUT(n,1)  - ( sumLogPI_OUT(n,1)  );
        PI_IN(n,t) = LAMBDA_FWD(n,1) * LAMBDA_BWD(n,1) / ( (1-LAMBDA_FWD(n,1))*(1-LAMBDA_BWD(n,1))*exp(tmp) +  LAMBDA_FWD(n,1) * LAMBDA_BWD(n,1) );
    end
end


sumProd_PI_OUT = prod(PI_OUT, 2);
sumProd_PI_OUT = max(min(sumProd_PI_OUT, 1-epsilon), epsilon);
sumProd1MinusPI_OUT = prod((1-PI_OUT), 2);
sumProd1MinusPI_OUT = max(min(sumProd1MinusPI_OUT, 1-epsilon), epsilon);

%PI_IN = NaN*ones(N,T);
% 
% % First the forward pass
% for n = 2 : 1 : N
%     LAMBDA_FWD(n,1) = ( P10(n,1) .* sumProd1MinusPI_OUT(n-1,1) .* (1 - LAMBDA_FWD(n-1,1)) ...
%         + (1 - P01(n,1)) .* sumProd_PI_OUT(n-1,1) .* LAMBDA_FWD(n-1,1) ) ./ ...
%     ( sumProd1MinusPI_OUT(n-1,1) .* (1 - LAMBDA_FWD(n-1,1)) ...
%     + sumProd_PI_OUT(n-1,1) .* LAMBDA_FWD(n-1,1) );
% end

% for n = 2 : 1 : N
%     LAMBDA_FWD(n,1) = ( P10(n,1) * (1 - LAMBDA_FWD(n-1,1)) * exp(sumLog1MinusPI_OUT(n-1,1))...
%         + (1 - P01(n,1)) * exp(sumLogPI_OUT(n-1,1)) * LAMBDA_FWD(n-1,1) ) / ...
%     ( (1 - LAMBDA_FWD(n-1,1)) * exp(sumLog1MinusPI_OUT(n-1,1)) ...
%     + exp(sumLogPI_OUT(n-1,1)) * LAMBDA_FWD(n-1,1) );
% end
% for n = 2 : 1 : N
%     LAMBDA_FWD(n,1) = ( P10(n,1) * (1 - LAMBDA_FWD(n-1,1)) ...
%         + (1 - P01(n,1)) * exp(sumLogPI_OUT(n-1,1)-sumLog1MinusPI_OUT(n-1,1)) * LAMBDA_FWD(n-1,1) ) / ...
%     ( (1 - LAMBDA_FWD(n-1,1))  ...
%     + exp(sumLogPI_OUT(n-1,1)-sumLog1MinusPI_OUT(n-1,1)) * LAMBDA_FWD(n-1,1) );
% end

% % Now the backward pass
% for n = N-1 : -1 : 1
%     LAMBDA_BWD(n,1) = ( P01(n,1) .* sumProd1MinusPI_OUT(n+1,1) .* (1 - LAMBDA_BWD(n+1,:)) ...
%         + (1 - P01(n,1)) .* sumProd_PI_OUT(n+1,1) .* LAMBDA_BWD(n+1,:) ) ./ ...
%         ((1 - P10(n,1) + P01(n,1)) .* sumProd1MinusPI_OUT(n+1,1) .* (1 - LAMBDA_BWD(n+1,1)) + ...
%         (1 - P01(n,1) + P10(n,1)) .* sumProd_PI_OUT(n+1,1) .* LAMBDA_BWD(n+1,1) );
% end

% PI_IN calculation
% for n = 1:N
%     for t = 1:T
%         PI_IN(n,t) = (LAMBDA_FWD(n,1) .* LAMBDA_BWD(n,1) .* sumProd_PI_OUT(n,1) ./ PI_OUT(n,t)  ) ./ ( (1 - LAMBDA_FWD(n,1)) .* ...
%      (1 - LAMBDA_BWD(n,1)) .* sumProd1MinusPI_OUT(n,1) ./ PI_OUT(n,t)  +  LAMBDA_FWD(n,1) .* LAMBDA_BWD(n,1) .* sumProd_PI_OUT(n,1) ./ PI_OUT(n,t) );
%     end
% end

% PI_IN = (LAMBDA_FWD .* LAMBDA_BWD) ./ ((1 - LAMBDA_FWD) .* ...
%     (1 - LAMBDA_BWD) + LAMBDA_FWD .* LAMBDA_BWD);

% First compute posterior means
MU_S = (PI_OUT .* repmat(LAMBDA_FWD,1,T) .* repmat(LAMBDA_BWD,1,T) ) ./ ((1 - PI_OUT) .* ...
    (1 - repmat(LAMBDA_FWD,1,T)) .* (1 - repmat(LAMBDA_BWD,1,T)) + PI_OUT .* repmat(LAMBDA_FWD,1,T) .* ...
    repmat(LAMBDA_BWD,1,T));

PS0S0 = (1 - P10(1:N-1,1)) .* ( (1 - LAMBDA_FWD(1:N-1,1)) .* ...
    sumProd1MinusPI_OUT(1:N-1,1) ) .* ( (1 - LAMBDA_BWD(2:N,1)) .* ...
    sumProd1MinusPI_OUT(2:N,1) );
PS1S0 = P10(1:N-1,1) .* ( (1 - LAMBDA_FWD(1:N-1,1)) .* ...
    sumProd1MinusPI_OUT(1:N-1,1) ) .* ( LAMBDA_BWD(2:N,1) .* ...
    sumProd_PI_OUT(2:N,1) );
PS0S1 = P01(1:N-1,1) .* ( LAMBDA_FWD(1:N-1,1) .* ...
    sumProd1MinusPI_OUT(1:N-1,1) ) .* ( (1 - LAMBDA_BWD(2:N,1)) .* ...
    sumProd1MinusPI_OUT(2:N,1) );
PS1S1 = (1 - P01(1:N-1,1)) .* ( LAMBDA_FWD(1:N-1,:) .* ...
    sumProd1MinusPI_OUT(1:N-1,1) ) .* ( LAMBDA_BWD(2:N,:) .* ...
    sumProd_PI_OUT(2:N,1) );

% PS0S0 = (1 - P10(1:N-1,:)) .* ((1 - repmat(LAMBDA_FWD(1:N-1,:),1,T) ) .* ...
%     (1 - PI_OUT(1:N-1,:))) .* ((1 - repmat(LAMBDA_BWD(2:N,:),1,T) ) .* ...
%     (1 - PI_OUT(2:N,:)));
% PS1S0 = P10(1:N-1,:) .* ((1 - repmat(LAMBDA_FWD(1:N-1,:),1,T) ) .* ...
%     (1 - PI_OUT(1:N-1,:))) .* ( repmat((LAMBDA_BWD(2:N,:)),1,T) .* ...
%     (PI_OUT(2:N,:)));
% PS0S1 = P01(1:N-1,:) .* (( repmat(LAMBDA_FWD(1:N-1,:),1,T) ) .* ...
%     (PI_OUT(1:N-1,:))) .* ((1 - repmat(LAMBDA_BWD(2:N,:),1,T) ) .* ...
%     (1 - PI_OUT(2:N,:)));
% PS1S1 = (1 - P01(1:N-1,:)) .* (( repmat(LAMBDA_FWD(1:N-1,:),1,T) ) .* ...
%     (PI_OUT(1:N-1,:))) .* (( repmat(LAMBDA_BWD(2:N,:),1,T) ) .* ...
%     (PI_OUT(2:N,:)));

S_CORR = PS1S1 ./ (PS0S0 + PS0S1 + PS1S0 + PS1S1);

% Now update the active-to-inactive transition probability via EM
P01_UPD = sum(sum(MU_S(1:N-1,:) - repmat(S_CORR,1,T))) / ...
    sum(sum(MU_S(1:N-1,:)));
% P01_UPD = sum(sum(MU_S(1:N-1,:) - S_CORR)) / ...
%     sum(sum(MU_S(1:N-1,:)));
P01_UPD = max(min(P01_UPD, 0.999), 0.001);
P01_UPD = repmat(P01_UPD,1,T);

S_POST = MU_S;


end