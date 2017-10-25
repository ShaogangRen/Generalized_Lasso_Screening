
clear all;
close all;

addpath('gls/');
addpath('data/');

fname = 'graph_sim_P3000_N_100_ER0.2_eps10';
load([fname '.mat']);

%%%%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
[N,P] = size(X);
mX = mean(X,2);
X = X - repmat(mX,[1,size(X,2)]);

gsolver = 'cvx';

[D, Node, FNode, Edge] = GraphNetwork(A);
uN = size(D,1);

T = TransformT(P, Edge);

tD = D*T;
tX = X*T;
ucell = get_ucell(tD,P);
tXCNm = sqrt(sum(tX.*tX,1))';

%%%%%%%%%%%%%%%% max lambda %%%%%%%%%%%%%%%%%%%
tXY = abs(tX'*Y);
[lamMax1, xImax ]= max(tXY);

XtY = X'*Y;
us1 = D*inv(D'*D)*XtY;  
lamMax2 = max(abs(us1));
lamMax = min([lamMax1, lamMax2]);


hlam = 0.001;
llam = 0.0005;
LamBuf = lamMax*[ 1 hlam : (-(hlam - llam)/50): llam]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Debug = true;
lamO =  LamBuf(1);  
betaO =  ones(P,1); 

BetaBuf = [betaO];
NLam = length(LamBuf);
RmRateBuf = [];

timebuf = [];
DSBuf = [];

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
node_list = [];
for i = 1:P
    nn = length(Node{i});
    ii = i*ones(nn,1);
    node_list = [node_list; ii, Node{i}'];
end

edge_list = [];
for i = 1:uN
    edge_list = [edge_list; Edge(i,:), 1, -1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('start the main algorithm!!');

GraphType = 0;
if uN == (P-1)
    GraphType = 1; %tree
elseif uN > (P-1)
    GraphType = 2; %graph with loop
elseif uN < (P-1)
    GraphType = 3; % forest   
end

if GraphType == 1
 Bound = 1; % start with good estimation
else 
 Bound = 2;
end

YNorm = norm(Y,2);

st = cputime;
for ilam = 2:NLam
       
        st1 = cputime;
        lam = LamBuf(ilam);
        betaO = BetaBuf(:, ilam-1);
        lamO =  LamBuf(ilam -1);

        %%%%%%%%%%%%%%% projection %%%%%%%%%%%%%%%%%%%%
        if Bound == 1   
            thetaO = Y/lamO;
            v = tX(:,xImax);
            v1 = sign(v'*thetaO)*v;
            
            v2 = Y/lam - thetaO;
            v1 = v1/norm(v1);
            Pv2 = v2 - v1*((v1'*v2));

            rr = norm(Pv2)/2;
            br = rr*tXCNm;
            cter = tX'*(thetaO + Pv2/2);
            L = cter - br;
            H = cter + br;
        end
        
        
        if Bound == 2
          thetaO = (Y - X*betaO)/lamO;
          lamDiff = (1/lam - 1/lamO);
          cy = lamDiff*Y/2;
          rr = abs(lamDiff)*YNorm/2;
          br = rr*tXCNm;
          cter = tX'*(thetaO + cy);
          L = cter - br;
          H = cter + br;
        end
        
        if Bound == 3
            thetaO = (Y - X*betaO)/lamO;
            v1 = Y/lamO - thetaO;
            v2 = Y/lam - thetaO;
            v1 = v1/norm(v1);
            Pv2 = v2 - v1*((v1'*v2));

            rr = norm(Pv2)/2;
            br = rr*tXCNm;
            cter = tX'*(thetaO + Pv2/2);
            L = cter - br;
            H = cter + br;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t1 = cputime - st1;
        st2 = cputime;
        
        u_p = ones(uN,1); %%u plus
        u_m = ones(uN,1); %%u minus
        [rmlab,u_p,u_m] = tBP_Edge(ucell,P,u_p,u_m,L,H);  %????

        reject_rate_i = sum(rmlab)/length(rmlab);
        RmRateBuf = [ RmRateBuf, reject_rate_i];
        t2 = cputime - st2;
       
        st3 = cputime;
        dbrmlab = double(rmlab);
        [newnode, newD]=NetShrink(dbrmlab, node_list, edge_list);
        newP = newnode(end,1); 
        newFNode = cell(newP,1); 
        
        rX = [];
        for i = 1:newP
            xids = newnode(newnode(:,1) == i,2);
            newFNode{i} = xids;
            rX = [rX, sum(X(:,xids),2)];
        end
        t3 = cputime - st3;
        
        st4 = cputime;
        if strcmp(gsolver, 'cvx')
            [tbeta,cvx_opt1] =  gfl_solver_cvx(rX,Y,lam,newD);
        end
           
        if Bound ~= 3
            if Bound == 1
                Bound = 3;
            end
            if Bound == 2
                if length(unique(tbeta)) ~= 1
                    Bound = 3;
                end
            end
        end
        
        beta = zeros(P,1);
        for i = 1:newP
           beta(newFNode{i}) = tbeta(i);
        end
        t4 = cputime - st4;
        
        BetaBuf = [BetaBuf, beta];
        timebuf = [timebuf; t1 t2 t3 t4];
        disp(['=====i=' num2str(ilam) '=== lam=' num2str(lam) '==== rej_rate=' ...
               num2str(reject_rate_i) '=========']);
end

usedtime = cputime - st
cvx_time = sum(timebuf(:,4))
lprogtime = sum(timebuf(:,2))
graphtime = sum(timebuf(:,3))
avg_rej_rate = mean(RmRateBuf)
%TBetaBuf = T*BetaBuf;
%sfname = [ 'gflscreen_' fname '.mat'];
%save(sfname,'LamBuf','BetaBuf','DSBuf', 'timebuf', 'RmRateBuf', 'usedtime' , 'cvx_time', 'lprogtime', 'graphtime', 'avg_rmrate');


