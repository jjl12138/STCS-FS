function [ xA_pri, vA_pri ] = B2A_Cal( xB_post, vB_post, xB_pri, vB_pri )
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
    vB2A_ext = 1/(1/vB_post - 1/vB_pri);
    vB2A_ext = min(1e8, max(1e-8, vB2A_ext));
    xB2A_ext = vB2A_ext/vB_post*xB_post - vB2A_ext/vB_pri*xB_pri;

    vA_pri = vB2A_ext;
    xA_pri = xB2A_ext;
end

