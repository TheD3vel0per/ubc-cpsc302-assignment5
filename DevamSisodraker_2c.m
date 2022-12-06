function J_GS(n)
close all
gcf
% Apply a stationary method to a linear system involving the n^2-by-n^2 
% Laplacian 

A=delsq(numgrid('S',n+2));
D=diag(diag(A));   % M=D for Jacovi
E=tril(A);         % M=E for Gauss-Seidel
itermax=10000;    % maximum number of iterations
resvecJ=zeros(itermax,1);  % allocate initial space for Jacobi residual norm vector 
resvecGS=zeros(itermax,1); % allocate initial space for Gauss-Seidel residual vector 
resvecSOR=zeros(itermax,1); % allocate initial space for SOR residual vector 

% Jacobi
xJ=zeros(n^2,1);
b=A*ones(n^2,1);   % generate a solution of all 1s and a right-hand side
nb=norm(b);
r=b;
for i=1:10000
    resvecJ(i)=norm(r)/nb;  % relative residual norm
    if resvecJ(i)<1e-6, break,end  % terminate loop if stopping criterion is satisfied
    xJ=xJ+D\r;              % next iterate
    r=b-A*xJ;               % update residual
end
semilogy(resvecJ);          % plot relative residual norm for Jacobi

% Gauss-Seidel
xGS=zeros(n^2,1);
r=b;
for i=1:10000
    resvecGS(i)=norm(r)/nb; % relative residual norm
    if resvecGS(i)<1e-6, break,end  % terminate loop if stopping criterion is satisfied
    xGS=xGS+E\r;            % next iterate
    r=b-A*xGS;              % update residual
end
hold on
semilogy(resvecGS,'r');     % plot relative residual norm for Gauss-Seidel

% SOR
xSOR=zeros(n^2,1);
r=b;
for i=1:10000
    wOpt = 2 / (1 + sin(pi / (n + 1)));
    resvecSOR(i)=norm(r)/nb; % relative residual norm
    if resvecSOR(i)<1e-6, break,end  % terminate loop if stopping criterion is satisfied
    xSOR=xSOR+((1 - wOpt)*D + wOpt*E)\r;            % next iterate
    r=b-A*xSOR;              % update residual    
end
hold on
semilogy(resvecSOR,'g');     % plot relative residual norm for Gauss-Seidel


legend('Jacobi','Gauss-Seidel','SOR')
saveas(gcf, "DevamSisodraker_2c.jpg", "jpg");