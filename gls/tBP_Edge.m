function [rmlab,u_p,u_m] = tBP_Edge(ucell,P,u_p,u_m,L,H)
 maxST = 10;

 for t = 1:maxST
 for i = 1:P
       if ~isempty(ucell{i})
                iUArr = ucell{i}(1,:);
                iUSgn = ucell{i}(2,:);
                iEdgeN = length(iUArr);

                lBuf = [ ]; % lasso penalty u_m(i)*LamRat
                hBuf = [ ]; % u_p(i)*LamRat

                for j = 1:iEdgeN    % graph penalty 
                    iu = iUArr(j);
                 %   ieu = ie + P;
                    if iUSgn(j) > 0
                         lBuf = [ lBuf u_m(iu)]; 
                         hBuf = [ hBuf u_p(iu)];
                    else
                         lBuf = [ lBuf u_p(iu)]; 
                         hBuf = [ hBuf u_m(iu)];
                    end
                end
                Sl = sum(lBuf) - L(i);
                Sh = sum(hBuf) + H(i);

                %ubuf = [(Sl - lBuf(1))/LamRat, (Sh - hBuf(1))/LamRat; u_p(i),
                %u_m(i)];[  
                
                ubuf = [];
                for j = 1:iEdgeN
                     iu = iUArr(j);
                     if iUSgn(j) > 0
                       % ieu = ie + P;
                        tbuf = [ Sl - lBuf(j), Sh - hBuf(j);   u_p(iu),  u_m(iu)];
                     else
                       % ie = -ie;  %ieu = ie + P;
                        tbuf = [ Sh - hBuf(j),  Sl - lBuf(j);  u_p(iu),  u_m(iu)];
                     end
                     ubuf = [ubuf tbuf ];
                end

                uval = min(ubuf);
                tvec = [ -ones(size(uval)); uval ];
                uval = max(tvec);
               % ieBuf = abs(ieBuf);

                for j = 1:iEdgeN
                  iu = iUArr(j);  %ieu = ie + P;
                  tt = (j-1)*2 + 1;
                  u_p(iu) = uval(tt);
                  u_m(iu) = uval(tt + 1);
                end
       end
 end
 end

%  trmN = (sum(u_p<1) + sum(u_m < 1));
 mg = 1 - 0.01;
 rmlab = (u_p < mg) &  (u_m < mg);





