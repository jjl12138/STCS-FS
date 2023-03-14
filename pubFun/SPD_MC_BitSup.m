function [ PI_IN, S_POST, P01_UPD,  LAMBDAT_UPD] = SPD_MC_BitSup( PI_OUT, LAMBDA, p01, LAMBDAT )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

%   forward and backward algorithm
% binary symbols
% -1 1
% 0  1
N = size(PI_OUT,1);
T = size(PI_OUT,2);

%P01 = p01*ones(N,T);
P01 = repmat(p01,N,1);
LAMBDA = repmat(LAMBDA,N,1);

P10 = P01 .* (LAMBDA ./ (1 - LAMBDA));

LAMBDA_FWD = NaN*ones(N,T);
LAMBDA_BWD = NaN*ones(N,T);

LAMBDA_FWD2 = NaN*ones(N,T);
LAMBDA_BWD2 = NaN*ones(N,T);

LAMBDA_FWD(1,:) = LAMBDA(1,:);
LAMBDA_BWD(N,:) = 0.5;

LAMBDA_T = LAMBDAT;

THETA_OUT = repmat(LAMBDA_T,N,1);
%THETA_IN = NaN*ones(N,T);
PT1 = NaN*ones(N,T);
PT0 = NaN*ones(N,T);

EPSILON = 1e-11;
Iter = 5;

for i = 1:Iter

% the forward pass for LAMBDA_FWD
% for n = 2 : 1 : N
%     LAMBDA_FWD(n,:) = (P10(n,:) .* (1 - PI_OUT(n-1,:)) .* (1 - LAMBDA_FWD(n-1,:)) + ...
%         (1 - P01(n,:)) .* PI_OUT(n-1,:) .* LAMBDA_FWD(n-1,:)) ./ ...
%         ((1 - PI_OUT(n-1,:)) .* (1 - LAMBDA_FWD(n-1,:)) + ...
%         PI_OUT(n-1,:) .* LAMBDA_FWD(n-1,:));
% end

% the forward pass for LAMBDA_FWD2, LAMBDA_FWD
LAMBDA_FWD(1,:) = THETA_OUT(1,:).*LAMBDA(1,:) + (1-THETA_OUT(1,:)).*EPSILON;
for n = 1 : 1 : N
    LAMBDA_FWD2(n,:) = PI_OUT(n,:).*LAMBDA_FWD(n,:)./( PI_OUT(n,:).*LAMBDA_FWD(n,:) ...
        + (1-PI_OUT(n,:)).*(1-LAMBDA_FWD(n,:)));
    
    if n<N
    LAMBDA_FWD(n+1,:) = (1 - P01(n,:)).*LAMBDA_FWD2(n,:).*THETA_OUT(n+1,:) + ...
        P10(n,:).*(1 - LAMBDA_FWD2(n,:)).*THETA_OUT(n+1,:) + EPSILON.*(1-THETA_OUT(n+1,:));
    end
    
end

% % Now the backward pass
% for n = N-1 : -1 : 1
%     LAMBDA_BWD(n,:) = (P01(n,:) .* (1 - PI_OUT(n+1,:)) .* (1 - LAMBDA_BWD(n+1,:)) + ...
%         (1 - P01(n,:)) .* PI_OUT(n+1,:) .* LAMBDA_BWD(n+1,:)) ./ ...
%         ((1 - P10(n,:) + P01(n,:)) .* (1 - PI_OUT(n+1,:)) .* (1 - LAMBDA_BWD(n+1,:)) + ...
%         (1 - P01(n,:) + P10(n,:)) .* PI_OUT(n+1,:) .* LAMBDA_BWD(n+1,:));
% end

% Now the backward pass for LAMBDA_BWD, LAMBDA_BWD2
%LAMBDA_BWD(N,:) = THETA_OUT(N+1,:).*(1/2) + (1-THETA_OUT(N+1,:)).*EPSILON;
for n = N : -1 : 1
    if n<N
    LAMBDA_BWD(n,:) = ( LAMBDA_BWD2(n+1,:).*THETA_OUT(n+1,:).*(1-P01(n,:)) ...
        + LAMBDA_BWD2(n+1,:).*(1 - THETA_OUT(n+1,:)).*EPSILON ...
        + (1-LAMBDA_BWD2(n+1,:)).*THETA_OUT(n+1,:).*P01(n,:) ...
        + (1-LAMBDA_BWD2(n+1,:)).*(1 - THETA_OUT(n+1,:)).*(1-EPSILON) )...
        ./( 2*EPSILON.*(1-THETA_OUT(n+1,:)).*LAMBDA_BWD2(n+1,:) ...
        + 2*(1-EPSILON).*(1-THETA_OUT(n+1,:)).*(1-LAMBDA_BWD2(n+1,:)) ...
        + LAMBDA_BWD2(n+1,:).*THETA_OUT(n+1,:).*(1+P10(n,:)-P01(n,:)) ...
        + (1-LAMBDA_BWD2(n+1,:)).*THETA_OUT(n+1,:).*(1+P01(n,:)-P10(n,:)) );
    end
    
    LAMBDA_BWD2(n,:) = PI_OUT(n,:).*LAMBDA_BWD(n,:)./( PI_OUT(n,:).*LAMBDA_BWD(n,:) ...
        + (1-PI_OUT(n,:)).*(1-LAMBDA_BWD(n,:)) );
end

% message passing for THETA_IN
PT1(1,:) = LAMBDA_BWD2(1,:).*LAMBDA(1,:) + (1-LAMBDA_BWD2(1,:)).*(1-LAMBDA(1,:));
PT0(1,:) = LAMBDA_BWD2(1,:).*EPSILON + (1-LAMBDA_BWD2(1,:)).*(1-EPSILON);
%PT1(N+1,:) = LAMBDA_FWD2(N,:).*(1/2) + (1-LAMBDA_FWD2(N,:)).*(1-1/2);
%PT0(N+1,:) = LAMBDA_FWD2(N,:).*EPSILON + (1-LAMBDA_FWD2(N,:)).*(1-EPSILON);
for n = 2 : 1 : N
    PT1(n,:) = LAMBDA_FWD2(n-1,:).*LAMBDA_BWD2(n,:).*(1-P01(n,:)) ...
        + LAMBDA_FWD2(n-1,:).*(1-LAMBDA_BWD2(n,:)).*P01(n,:) ...
        + (1-LAMBDA_FWD2(n-1,:)).*LAMBDA_BWD2(n,:).*P10(n,:) ...
        + (1-LAMBDA_FWD2(n-1,:)).*(1-LAMBDA_BWD2(n,:)).*(1-P10(n,:));
    PT0(n,:) = (1-LAMBDA_BWD2(n,:)).*(1-EPSILON) + LAMBDA_BWD2(n,:).*EPSILON;
end
THETA_IN = PT1./(PT1+PT0);

% message passing for THETA_OUT
% SUM_THETA_IN = prod(THETA_IN,1);
% SUM_1MINUS_THETA_IN = prod(1-THETA_IN,1);
% 
% for n = 1 : 1 : N+1
%     THETA_OUT(n,:) = ( LAMBDA_T.*SUM_THETA_IN./THETA_IN(n,:) )./( LAMBDA_T.*SUM_THETA_IN./THETA_IN(n,:) ...
%         + (1-LAMBDA_T).*SUM_1MINUS_THETA_IN./(1 - THETA_IN(n,:)) );
% end
SUM_LLR_THETA_IN = sum(log((1-THETA_IN)./THETA_IN),1);
LLR_LAMBDA_T = log((1-LAMBDA_T)./LAMBDA_T);
for n = 1 : 1 : N
    THETA_OUT(n,:) = 1./(1 + exp(LLR_LAMBDA_T + SUM_LLR_THETA_IN - log( (1-THETA_IN(n,:))./THETA_IN(n,:) )));
end


end
% LOOP Iter end

% Update LMABDAT via EM 
LAMBDAT_UPD = 1./(1 + exp(LLR_LAMBDA_T + SUM_LLR_THETA_IN)); 

% message passing for PI_IN
PI_IN = (LAMBDA_FWD .* LAMBDA_BWD) ./ ((1 - LAMBDA_FWD) .* ...
    (1 - LAMBDA_BWD) + LAMBDA_FWD .* LAMBDA_BWD);

% First compute posterior means
MU_S = (PI_OUT .* LAMBDA_FWD .* LAMBDA_BWD) ./ ((1 - PI_OUT) .* ...
    (1 - LAMBDA_FWD) .* (1 - LAMBDA_BWD) + PI_OUT .* LAMBDA_FWD .* ...
    LAMBDA_BWD);

PS0S0 = (1 - P10(1:N-1,:)) .* ((1 - LAMBDA_FWD(1:N-1,:)) .* ...
    (1 - PI_OUT(1:N-1,:))) .* ((1 - LAMBDA_BWD(2:N,:)) .* ...
    (1 - PI_OUT(2:N,:)));
PS1S0 = P10(1:N-1,:) .* ((1 - LAMBDA_FWD(1:N-1,:)) .* ...
    (1 - PI_OUT(1:N-1,:))) .* ((LAMBDA_BWD(2:N,:)) .* ...
    (PI_OUT(2:N,:)));
PS0S1 = P01(1:N-1,:) .* ((LAMBDA_FWD(1:N-1,:)) .* ...
    (PI_OUT(1:N-1,:))) .* ((1 - LAMBDA_BWD(2:N,:)) .* ...
    (1 - PI_OUT(2:N,:)));
PS1S1 = (1 - P01(1:N-1,:)) .* ((LAMBDA_FWD(1:N-1,:)) .* ...
    (PI_OUT(1:N-1,:))) .* ((LAMBDA_BWD(2:N,:)) .* ...
    (PI_OUT(2:N,:)));

S_CORR = PS1S1 ./ (PS0S0 + PS0S1 + PS1S0 + PS1S1);

% Now update the active-to-inactive transition probability via EM
P01_UPD = NaN*ones(1,T);
for t = 1:T
    P01_UPD(1,t) = sum(sum(MU_S(1:N-1,t) - S_CORR(:,t))) / ...
        sum(sum(MU_S(1:N-1,t)));
end

P01_UPD = max(min(P01_UPD, 1), 0);

S_POST = MU_S;


end

