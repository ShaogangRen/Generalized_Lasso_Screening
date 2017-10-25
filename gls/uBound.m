function [rmlab,uL, uH] = uBound(ucell,P,uL,uH,L,H)

 maxST = 10;
 
 DValP = cell(P,1);
 DValN = cell(P,1);
 DValI = cell(P,1);
 
for i = 1:P
       if ~isempty(ucell{i})
           
                iUArr = ucell{i}(1,:);
                iDVal = ucell{i}(2,:);
                %iEdgeN = length(iUArr);
                
                pidx = (iDVal >0);
                DValI{i} = pidx;
                vbuf = zeros(size(iDVal));
                vbuf(pidx) = iDVal(pidx);
                DValP{i} = vbuf;
                vbuf = zeros(size(iDVal));
                nidx = ~pidx;
                vbuf(nidx) = iDVal(nidx);
                DValN{i} = vbuf;
       end 
end
 


 for t = 1:maxST
      flag = 0;
 for i = 1:P
    
     
       if ~isempty(ucell{i})
           
                iUArr = ucell{i}(1,:);
                iDVal = ucell{i}(2,:);
                iEdgeN = length(iUArr);
                
                
                MaBuf = DValP{i}.*uH(iUArr)' + DValN{i}.*uL(iUArr)';
                MiBuf = DValP{i}.*uL(iUArr)' + DValN{i}.*uH(iUArr)';
                
                duMax = sum(MaBuf);
                duMin = sum(MiBuf);
                pidx = DValI{i};
                
                
                for j = 1:iEdgeN
                    
                    iu = iUArr(j);
                    if pidx(j) 
                        uht = (H(i) - duMin + MiBuf(j))/iDVal(j);
                        ult = (L(i) - duMax  + MaBuf(j))/iDVal(j);
                    else
                        ult = (H(i) - duMin + MiBuf(j))/iDVal(j);
                        uht = (L(i) - duMax  + MaBuf(j))/iDVal(j);
                    end
                    
                    if uht < uH(iu)
                        uH(iu) = uht;
                        flag = 1;
                    end
                    
                                        
                    if ult > uL(iu)
                        uL(iu) = ult;
                        flag = 1;
                    end
                    
                end

       end
       
       
 end
 
     if flag == 0
       break; 
     end
 end

 mg = 1 - 0.0001;
 rmlab = (uH < mg) & (uL > (-mg));




