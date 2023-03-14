function [ x, xt, rp ] = MultiChanGen1( N, P, LTimDom, KTimDom, KAngDom, grpNumAngDom )
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明
    
%     N = 256;
%     P = 64;
%     KAngDom = 64;
%     grpNumAngDom = 4;
%     LTimDom = 32;
%     KTimDom = 4;

    xtmptim = zeros(N,P);
    xtmpfre = zeros(N,P);
    x = cell(1,P);
    xt = cell(1,P);
    
    Nsmal = N/grpNumAngDom;
    grpsmal = KAngDom/grpNumAngDom;
    
    xgrptmp = zeros(N,grpNumAngDom);
    for ii = 1:grpNumAngDom
        %loc = randi([1 Nsmal-grpsmal+1],1,1);
        loc = randi([1 Nsmal-grpsmal],1,1);
        loc = loc + (ii-1)*Nsmal;
        xgrptmp(loc:loc+grpsmal-1,ii) = ones(grpsmal,1);
        
    end
    
    rp = randperm(LTimDom);
    rp = rp(1:KTimDom);
    for ii = 1:KTimDom
        iii = mod(ii-1,grpNumAngDom)+1;
        index = xgrptmp(:,iii)==1;
        xtmptim(index, rp(1,ii)) = 1; 
    end
    
    xtmptim = xtmptim.*(randn(N,P)+1i*randn(N,P))/sqrt(2);
    xtmptim = xtmptim/sqrt(var(xtmptim(:)));
    
    F = dftmtx(P)/sqrt(P);
    xtmpfre = xtmptim * F.';
%     for ii = 1:KAngDom
%         xtmpfre(index(ii,1), :) = xtmptim(index(ii,1), :) * F.';
%     end
    %norScal = sqrt(var(xtmpfre(:)));
    %x = xtmpfre/norScal;
    norScal = 1;
    
    for ii = 1:P
        x{ii} = xtmpfre(:,ii)/norScal;
    end
    
%    xt1 = xtmpfre/norScal*(F').';
%     for ii = 1:P
%         xt{ii} = xt1(:,ii);
%     end
    for ii = 1:P
        xt{ii} = xtmptim(:,ii);
    end

end