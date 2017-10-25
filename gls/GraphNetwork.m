function [D, Node, FNode, Edge] = GraphNetwork(A)
%dsize = [3,3];
%A = GraphAdj(dsize);
[I,J,V] = find(A);
D = [];
nN = size(A,1);
eN = length(I);

for i = 1:length(I)
    td = zeros(1,nN);
    td(I(i)) = 1;
    td(J(i)) = -1;
    ts = sparse(td);
%     size(D)
%     size(ts)
    D = [D; ts ];
end

D = sparse(D);

%%% generate network %%%
Node = cell(nN,1);
FNode = cell(nN,1);
Edge = zeros(eN,2);

for i = 1:eN
     Edge(i,:) = [I(i),J(i)]; 
     Node{I(i)} = [ Node{I(i)}, i];
     Node{J(i)} = [ Node{J(i)}, i];
end

for i = 1:nN
    FNode{i} = i;
end



