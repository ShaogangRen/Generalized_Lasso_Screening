function [rmlab, u_l, u_h] =  LP_BoundEst(uN, D, u_lb, u_hb, L, H)

u_h = zeros(uN,1);
u_l = zeros(uN,1);

dD = [D'; -D'];
b = [H; -L];

zf = zeros(uN,1);
for i = 1:uN
    tf = zf;
    tf(i) = 1;
    [uu, fu ] = linprog(tf,dD,b,[],[], u_lb, u_hb);
     u_l(i) = fu;
     
    tf(i) = -1;
    [uu, fu ] = linprog(tf,dD,b,[],[], u_lb, u_hb);
    u_h(i) = -fu;
     
end

 mg = 1 - 0.00001;
 rmlab = (u_l > -mg) &  (u_h < mg);




