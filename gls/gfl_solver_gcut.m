function [beta, funcVal] = gfl_solver_gcut(X,Y,rho1,rho2,D)
%A = A + A';

EN = size(D,1);
I = zeros(EN,1);
J = zeros(EN,1);

%[M,N,V] = find(D);

for i = 1:EN
%    if I(M(i)) ==0
%       I(M(i)) = N(i);
%    else
%       J(M(i)) = N(i);
%    end
    M = find(D(i,:));
    I(i) = M(1);
    J(i) = M(2);
end

Graph{1} = EN;                  
Graph{2} = ones(EN,1);   %weight         
Graph{3} = I;                  
Graph{4} = J;

opts.tol = 10^-6;   % tolerance.
opts.maxIter = 10000; % maximum iteration number of optimization.

tic
[beta, funcVal] = fast_gfl_rho(X, Y, Graph, rho1, rho2, opts);
t = toc;





