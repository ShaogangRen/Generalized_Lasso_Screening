function W = GraphAdj(dsize, eRate)
% numSites = prod(dsize);
% id1 = [1:numSites, 1:numSites];
% %id1 = [1:numSites, 1:numSites, 1:numSites];
% id2 = [ 1+1:numSites+1,...
%         1+dsize(1):numSites+dsize(1)]; %...
%         %1+sizeData(1)*sizeData(2):numSites+sizeData(1)*sizeData(2)];
% value = ones(1,2*numSites);
% W = sparse(id1,id2,value);
% W = W(1:numSites,1:numSites);

numSites = dsize(1)*dsize(2);
id1 = [1:numSites, 1:numSites];
id2 = [(1+1):(numSites+1),  ...
       1+ dsize(1) : numSites+dsize(1) ];
   
value = ones(1,2*numSites);

%%%%% remove invalid links
%for i = 1:dsize(3)
%     for j = 1: dsize(2)
%        t = j*dsize(1);
%        value(t) = 0;
%     end
%end
value((end - dsize(1) + 1) : end) = 0;

%remove edges
nMaxEdge = numSites - dsize(1);
%eRate = 0.2;
%EdgeIdx = numSites + randi(nMaxEdge,[1,round((1 - eRate)*nMaxEdge)]);
%EdgeIdx = numSites + randi(nMaxEdge,[1,round((1 - eRate)*nMaxEdge)]);
%value(EdgeIdx) = 0;

st = (numSites + numSites/20);
ed = st + numSites*(1-eRate);
value(st:ed) = 0;




W = sparse(id1,id2,value);
W = W(1:numSites,1:numSites);









