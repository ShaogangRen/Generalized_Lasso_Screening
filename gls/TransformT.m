function T = TransformT(P, edge_list)


[NodeData, PathNode]=GenTransfT(P, edge_list);

%NDBox = cell(P,1);
T = zeros(P,P);

for i = 1:P
     T(NodeData((NodeData(:,1) == i),2),i)=1;
end

   pT = [];
if size(PathNode,1) > 0
    PathN = PathNode(end,1);

    for i = 1:PathN
        path = PathNode((PathNode(:,1) == i), 2);
        tt = path(1);   flag = false;
        if tt == path(end)
           path = path(1:end-1); 
           flag = true;
        end
      
        plth = length(path);
        for j = 2: plth
            for k = 1:(j-1)
                tt = sum(T(:,path(k:j)),2);
                pT = [pT, tt];
            end  
        end
        
        if flag
           path = [ path(2:end); path(1) ];
           for k= 2:(plth-1)
              tt = sum(T(:,path(k:plth)),2);
              pT = [pT, tt];  
           end
        end
        
    end
end

T = [T, pT];



