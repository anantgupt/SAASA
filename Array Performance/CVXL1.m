function [x,info] = CVXL1(A,b,sigma,solver)
% Solve using CVX:
%
% Minimize ||x||_1 subject to ||Ax-b||_2 <= sigma
%

[m,n] = size(A);
info  = struct();

% Set solver
solver_old = cvx_solver;
if ~isempty(solver), cvx_solver(solver); end

% Start timer
tic;

if isreal(A) && isreal(B)
   cvx_begin
      variable x(n);
      minimize norm(x,1);
      subject to
         norm(A*x-b,2) <= sigma;
   cvx_end
else
   cvx_begin quiet
   cvx_precision low
      variable x(n) complex;
      minimize norm(x,1);
      subject to
         norm(A*x-b,2) <= sigma;
   cvx_end
end

% Stop timer
info.time_solve = toc;
info.time_prepare = 0;

% Restore solver
cvx_solver(solver_old);
