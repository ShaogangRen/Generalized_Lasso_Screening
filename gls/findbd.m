function [lowbd, highbd] = findbd(M, lbd, hbd)

[rr, cc] = size(M);
ps = M >0;
psM= zeros(rr,cc);
psM(ps) = M(ps);
%clear ps;
ngM =  zeros(rr,cc);
ngM(~ps) = M(~ps);
clear ps;

%hB = repmat(hbd,[rr,1]);
%lB = repmat(lbd,[rr,1]);

%BUF = zeros([rr, cc]);
%BUF(ps) = hB(ps);
%BUF(ng) = lB(ng);
highbd = psM*hbd + ngM*lbd;%sum(BUF .* M,2);

%BUF = zeros([rr, cc]);
%BUF(ps) = lB(ps);
%BUF(ng) = hB(ng);
lowbd = ngM*hbd + psM*lbd; % sum(BUF .* M,2);
