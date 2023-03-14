function [ estFin, optFin, estHist ] = tcsEst_mmv( y, A, opt, sigEst )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
estHist = [];
optFin = [];

nit = opt.nit;
%step = opt.step;
nuw = opt.nuw;
saveHist = opt.saveHist;
timSupp = opt.timSupp;
xtrue = opt.x; 
%tol = opt.tol;

T = size(A,2);
[M, N] = size(A{1});

xB_postPre = 0;

xA_pri = zeros(N,T);
vA_pri = ones(1,T);

if ~isempty(opt.xA_pri0)
    xA_pri = opt.xA_pri0;
end

if ~isempty(opt.vA_pri0)
    vA_pri = opt.vA_pri0;
end

if saveHist
    estHist.xA_post = nan(N*T,1);
    estHist.vA_post = nan(T,1);
    estHist.xB_pri = nan(N*T,1);
    estHist.vB_pri = nan(T,1);
    estHist.xB_post = nan(N*T,1);
    estHist.vB_post = nan(T,1);
    estHist.xA_pri = nan(N*T,1);
    estHist.vA_pri = nan(T,1);
    estHist.nuw = nan(N,1);
    estHist.errX = nan(1,1);
    estHist.difX = nan(1,1);
    estHist.difX2 = nan(1,1);
end

%gamwHat = 1/nuw;
% [U,D] = eig(A*A');
% d = diag(D);
% U = opt.U;
% D = opt.D;
% d = opt.d;
xB_pri = zeros(N,T);
vB_pri = zeros(1,T);

stop = false;
it = 0;
while ~stop
    
    % Iteration count
    it = it + 1;
    
    % Check for final iteration
    if it >= nit
        stop = true;
    end
    
    for t = 1:T
        [ xB_pri(:,t), vB_pri(:,t) ] = A2B_Cal( y{t}, A{t}, nuw{t}, xA_pri(:,t), vA_pri(:,t) );
    end
    
    [ xB_post, vB_post, ~] = CBG_Est( xB_pri, repmat(vB_pri, N, 1), sigEst);
    vB_post = mean(vB_post,1);
    
    if timSupp
        F = dftmtx(T)/sqrt(T);
        xBt_post = xB_post*(F').';
   
        % Soft Thresholding Denoiser
        [ xBt_post, ~ ] = CsoftThres_Est( xBt_post, sqrt(mean(vB_post,2)), []);
        % Cut-off Denoiser
        %xBt_post(:,LTimDom+1:1:end) = 0;

        xB_post = xBt_post*F.';

        [ xA_pri ] = Extrinsic_Est( xB_post, xB_pri, vB_pri, sigEst );

        for t = 1:T
            vA_pri(1,t) = ( norm(y{t}-A{t}*xA_pri(:,t))^2 - M*nuw{t})/M;
            vA_pri(1,t) = min(1e8, max(1e-8, vA_pri(1,t)));
        end       
    else
        for t = 1:T
            [ xA_pri(:,t), vA_pri(1,t) ] = B2A_Cal( xB_post(:,t), vB_post(1,t), xB_pri(:,t), vB_pri(1,t) );
        end 
    end
   
    
%     gam1overD = (1./d)*(1/vA_pri);
%     UyAr1 = U'*y - U'*A*xA_pri;
%     gam1overD_UyAr1_Sq = (gam1overD.^2).*abs(UyAr1).^2;
%     for ii=1:5
%     gamwHat_old = gamwHat;
%     % note that resNormSq = sum(abs(y-A*x1).^2,1);
%     resNormSq = sum(gam1overD_UyAr1_Sq./((gamwHat+gam1overD).^2)); 
%     gamwHat = M/mean(resNormSq + sum(1./(gamwHat+gam1overD),1) );
%     if abs(gamwHat_old-gamwHat)/gamwHat < 0.01, break; end;
%     end      
%     nuw = 1/gamwHat;
    
    % Check for convergence
    if (it>1) && (stop==false)
%         if ( norm(xB_postPre(:) - xB_post(:)) / norm(xB_post(:)) < tol )
%             stop = true;
%         end
    end
    
    
    if saveHist
        %estHist.xA_post(:,it) = xA_post(:);
        %estHist.vA_post(:,it) = vA_post(:);
        estHist.xA_post = zeros(N*T,1);
        estHist.vA_post = zeros(T,1);
        estHist.xB_pri(:,it) = xB_pri(:);
        estHist.vB_pri(:,it) = vB_pri(:);
        estHist.xB_post(:,it) = xB_post(:);
        estHist.vB_post(:,it) = vB_post(:);
        estHist.xA_pri(:,it) = xA_pri(:);
        estHist.vA_pri(:,it) = vA_pri(:);
        %estHist.nuw(:,it) = nuw;
        estHist.nuw(:,it) = zeros(N,1);
        estHist.errX(:,it) = (norm(xB_post-[xtrue{:}],'fro')/norm([xtrue{:}],'fro'))^2;
        estHist.difX(:,it) = norm(xB_postPre(:) - xB_post(:)) / norm(xB_post(:));
    end
    
    xB_postPre = xB_post;
    
end

estFin.xB_post = xB_post;
estFin.xB_pri = xB_pri;
estFin.vB_pri = vB_pri;
estFin.nuw = nuw;

estFin.xA_pri = xA_pri;
estFin.vA_pri = vA_pri;
estFin.it = it;


%optFin = opt;

end