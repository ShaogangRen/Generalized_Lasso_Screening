function [ newD, newFNode, newP, newEN, newNode, newEdge] = graph_reduce(rmlab, P, Node, FNode,Edge)
%% network shrinkage algorithm
%%%% rmlab, binary vector : 1 remove, 0 not remove

        rmeA = rmlab;
        EN = length(rmlab);
         %%% remove fused nodes
          rmvC  = false(1,P);
%          trmv = rmvA | rmvB;
         for ie = 1:EN
              if  rmeA(ie)
                 tEdge = Edge(ie,:);
                  eh = tEdge(1); et = tEdge(2); 
                % if ~(trmv(eh) |  trmv(et))
                     %% merge the data, then link neibors of eh to et
                     FNode{et} = [ FNode{et},  FNode{eh}];
                     rNode = Node{eh};
                     rlEN = length(rNode);

                      tm = Edge(Node{et},:)';  %matlab takes column as vector
                      tv = tm((tm ~= et));
                      tv = tv(tv ~= eh);
                      tv = [tv; et];
                      
                      tm = Edge(rNode,:)';
                      vc =  tm((tm ~= eh));
                     for irle = 1:rlEN
                         if all(~(vc(irle) == tv))
                            Edge(rNode(irle),Edge(rNode(irle),:) == eh) = et;
                            Node{et} = [Node{et}, rNode(irle)];
                         end
                     end
                     Node{eh} = []; rmvC(eh) = 1;
                     Node{et}(Node{et} == ie) = [];
                % else
                %     disp([ 'Debug: eh=' num2str(eh), ' et=' num2str(et)]);
                % end
              end
         end
         
         %%% remove nodes
          rmv = rmvC;
          newP = P - sum(rmv);
          newLabel = zeros(P,1);
          newLabel(~rmv) = [1:newP];
          newNode = cell(newP,1);
          newEdge =  [];
          newEN = 0;
          newD = [];
          
          for i = 1:P
              if newLabel(i) > 0
                  iNode = Node{i}; en = length(iNode);
                  for j = 1:en
                     tEdge = Edge(iNode(j),:);
                     if tEdge(1) == i
                         tnd = tEdge(2);
                         if newLabel(tnd)> 0
                            newEN = newEN + 1;
                            if  newLabel(i) < newLabel(tnd)
                                 ehead = newLabel(i); etail = newLabel(tnd);
                            else
                                 ehead = newLabel(tnd); etail = newLabel(i);
                            end

                            newEdge = [newEdge; ehead, etail];
                            newNode{newLabel(i)} = [newNode{newLabel(i)}, newEN];
                            newNode{newLabel(tnd)} = [newNode{newLabel(tnd)}, newEN];
                            td = sparse(1,newP);  td(ehead) = 1;  td(etail) = -1;
                            newD = [newD; td ];
                         end
                     end
                  end
              end
          end
        %  clear Node;
        %  clear Edge;
        FNode(rmv) = [];
        %  Node = newNode;
        %  Edge = newEdge;
        newFNode = FNode;
          
          
          
