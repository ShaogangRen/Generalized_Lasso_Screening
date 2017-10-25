function [rmlab,u_p,u_m] = BP_Edge(P,Node,Edge,u_p,u_m,L,H)

  for i = 1:P
        iNode = Node{i};
        iEdgeN = length(iNode);
        lBuf = [ ]; % lasso penalty u_m(i)*LamRat
        hBuf = [ ]; % u_p(i)*LamRat

        ieBuf = [];
        for j = 1:iEdgeN    % graph penalty 
            ie = iNode(j);
         %   ieu = ie + P;
            if i == Edge(ie,1)
                 lBuf = [ lBuf u_m(ie)]; 
                 hBuf = [ hBuf u_p(ie)];
            else
                 lBuf = [ lBuf u_p(ie)]; 
                 hBuf = [ hBuf u_m(ie)];
                 ie = -ie;
            end
            ieBuf = [ieBuf,ie];
        end
        Sl = sum(lBuf) - L(i);
        Sh = sum(hBuf) + H(i);

        %ubuf = [(Sl - lBuf(1))/LamRat, (Sh - hBuf(1))/LamRat; u_p(i),
        %u_m(i)];[
        ubuf = [];
        for j = 1:iEdgeN
             ie = ieBuf(j);  
             if ie > 0
               % ieu = ie + P;
                tbuf = [ Sl - lBuf(j), Sh - hBuf(j);   u_p(ie),  u_m(ie)];
             else
                ie = -ie;  %ieu = ie + P;
                tbuf = [ Sh - hBuf(j),  Sl - lBuf(j);  u_p(ie),  u_m(ie)];
             end
             ubuf = [ubuf tbuf ];
        end

        uval = min(ubuf);
        ieBuf = abs(ieBuf);

        for j = 1:iEdgeN
          ie = ieBuf(j);  %ieu = ie + P;
          tt = (j-1)*2 + 1;
          u_p(ie) = uval(tt);
          u_m(ie) = uval(tt + 1);
        end
  end

%  trmN = (sum(u_p<1) + sum(u_m < 1));
 rmlab = (u_p < 1) &  (u_m < 1);


