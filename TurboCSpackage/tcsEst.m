function [ estFin, optFin, estHist ] = tcsEst( y, A, opt, sigEst )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
estHist = [];
optFin = [];

nit = opt.nit;
step = opt.step;
nuw = opt.nuw;
saveHist = opt.saveHist;
xtrue = opt.x; 
tol = opt.tol;

[M, N] = size(A);

vA2BPre = 0;
vB2APre = 0;
xA2BPre = 0;
xB2APre = 0;
xB_postPre = 0;

xA_pri = zeros(N,1);
vA_pri = 1;

if ~isempty(opt.xA_pri0)
    xA_pri = opt.xA_pri0;
end

if ~isempty(opt.vA_pri0)
    vA_pri = opt.vA_pri0;
end

if saveHist
    estHist.xA_post = nan(N,1);
    estHist.vA_post = nan(N,1);
    estHist.xB_pri = nan(N,1);
    estHist.vB_pri = nan(N,1);
    estHist.xB_post = nan(N,1);
    estHist.vB_post = nan(N,1);
    estHist.xA_pri = nan(N,1);
    estHist.vA_pri = nan(N,1);
    estHist.nuw = nan(N,1);
    estHist.errX = nan(1,1);
    estHist.difX = nan(1,1);
    estHist.difX2 = nan(1,1);
end

gamwHat = 1/nuw;
% [U,D] = eig(A*A');
% d = diag(D);
U = opt.U;
D = opt.D;
d = opt.d;

stop = false;
it = 0;
while ~stop
    
    % Iteration count
    it = it + 1;
    
    % Check for final iteration
    if it >= nit
        stop = true;
    end
    
    xA_post = xA_pri + vA_pri / ( vA_pri + nuw ) * A' * ( y - A * xA_pri);
    vA_post = vA_pri - M/N * vA_pri^2 / (vA_pri + nuw);
    
    %update extrinsic
    %vA2B_ext = N/M * ( vA_pri + nuw) - vA_pri;
    vA2B_ext = 1/( 1/vA_post - 1/vA_pri );
    vA2B_ext = step*vA2B_ext + (1-step)*vA2BPre;

    vA2B_ext = min(1e11, max(1e-11, vA2B_ext));
    vA2BPre = vA2B_ext;
    
    xA2B_ext = vA2B_ext/vA_post*xA_post - vA2B_ext/vA_pri*xA_pri;
    xA2B_ext = step*xA2B_ext + (1-step)*xA2BPre;
    xA2BPre = xA2B_ext;
    
    xB_pri = xA2B_ext;
    vB_pri = vA2B_ext;
    
    [ xB_post, vB_post, ~] = CBG_Est( xB_pri, vB_pri, sigEst);
    vB_post = mean(vB_post);
    
    vB2A_ext = 1/(1/vB_post - 1/vB_pri);
    vB2A_ext = step*vB2A_ext + (1-step)*vB2APre;
    vB2A_ext = min(1e11, max(1e-11, vB2A_ext));
    vB2APre = vB2A_ext;
    
    xB2A_ext = vB2A_ext/vB_post*xB_post - vB2A_ext/vB_pri*xB_pri;
    

    %cOpt1 = real(xB_pri'* (xB_post - cOpt2*xB_pri) )/( (xB_post-cOpt2*xB_pri)'*(xB_post-cOpt2*xB_pri));
    %xB2A_ext = cOpt1*(xB_post - cOpt2*xB_pri);
%     fprintf('NMSE = %.4f\n', (norm(x - xtmp,'fro')/norm(x,'fro'))^2);
%     NMSE_bp_mea_tmp(it,nsim,snri) = (norm(x - xtmp,'fro')/norm(x,'fro'))^2;
%     xA_pri = xtmp;
    
    xB2A_ext = step*xB2A_ext + (1-step)*xB2APre;
%     if saveHist
%         estHist.difX2(:,it) = norm(xB2APre(:) - xB2A_ext(:)) / norm(xB2A_ext(:));
%     end
    xB2APre = xB2A_ext;
    
    %vB2A_ext = (norm(y-A*xB2A_ext)^2 - M*nuw)/M;
    vB2A_ext = step*vB2A_ext + (1-step)*vB2APre;
    vB2A_ext = min(1e11, max(1e-11, vB2A_ext));
    vB2APre = vB2A_ext;
    
    vA_pri = vB2A_ext;
    xA_pri = xB2A_ext;
    
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
        if ( norm(xB_postPre(:) - xB_post(:)) / norm(xB_post(:)) < tol )
            stop = true;
        end
    end
    
    
    if saveHist
        estHist.xA_post(:,it) = xA_post(:);
        estHist.vA_post(:,it) = vA_post(:);
        estHist.xB_pri(:,it) = xB_pri(:);
        estHist.vB_pri(:,it) = vB_pri(:);
        estHist.xB_post(:,it) = xB_post(:);
        estHist.vB_post(:,it) = vB_post(:);
        estHist.xA_pri(:,it) = xA_pri(:);
        estHist.vA_pri(:,it) = vA_pri(:);
        estHist.nuw(:,it) = nuw;
        estHist.errX(:,it) = (norm(xB_post-xtrue,'fro')/norm(xtrue,'fro'))^2;
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

