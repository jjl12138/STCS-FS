function [ xB_pri, vB_pri ] = A2B_Cal( y, A, nuw, xA_pri, vA_pri )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
    [M, N] = size(A);
    xA_post = xA_pri + vA_pri / ( vA_pri + nuw ) * A' * ( y - A * xA_pri);
    vA_post = vA_pri - M/N * vA_pri^2 / (vA_pri + nuw);

    %update extrinsic
    %vA2B_ext = N/M * ( vA_pri + nuw) - vA_pri;
    vA2B_ext = 1/( 1/vA_post - 1/vA_pri );
    vA2B_ext = min(1e8, max(1e-8, vA2B_ext));
    xA2B_ext = vA2B_ext/vA_post*xA_post - vA2B_ext/vA_pri*xA_pri;

    xB_pri = xA2B_ext;
    vB_pri = vA2B_ext;
end

