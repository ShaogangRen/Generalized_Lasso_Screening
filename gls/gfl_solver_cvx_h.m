function [beta, cvx_optval] = gfl_solver_cvx_h(X,Y,lam,D)

P = size(X,2);

cvx_begin quiet
    %cvx_precision high
    cvx_precision(3E-16)
    variable B(P,1)
    minimize(  0.5*(Y - X*B)'*(Y - X*B) + lam*norm(D*B,1) ) %0.5*square(
cvx_end

beta = B;

