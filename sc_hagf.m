function [H, C, Z] = sc_hagf(X, alpha, beta, eta)

% ||H - FX||_F^2 + alpha * ||H - CH||_F^2 + beta||C - FC||_F^2
% s.t. Ce = e, C = C', C>0, F = eta*F_l + (1- eta)*F_h, diag(C) = 0
% F_l = 3/4*I + 1/4*C, F_h = 1/4*I - 1/4*C

% X is a data matrix with size n x d, n is the number of sample, d is the
% original dimensionality, C is the reconstruction coefficient matrix





%% parameters
tol = 1e-7;
maxIter = 1200;
rho = 1.1;
max_mu = 1e30;
mu = 1e-6;

[n,~] = size(X);
e = ones(n,1);
I = eye(n);

XXt = X*X' ;
eet = e*e';

%% auxiliary Variables
H = X;
Z = zeros(n, n);


%% Lagrange multipliers
Y1 = zeros(n, n);
Y2 = zeros(n, 1);

%% Start main loop
iter  = 0;
while iter < maxIter
    iter = iter + 1;

    % update C
    HHt = H*H';
    I_Z = (3 - 2*eta)*I - (2*eta - 1)*Z;
    A = 2*(2*eta -1)^2*XXt + 2*alpha*HHt + 2*beta*(I_Z*I_Z') + mu*(I + eet);
    B = (2*eta - 1)*8*H*X' - 2*(4*eta^2 - 1)*XXt + 2*alpha*(HHt) + mu*(Z + eet) - Y1 - Y2*e';
    C = B/A;

    % update H
    I_C = I - C;
    A = alpha*(I_C'*I_C) + 16*I;
    B = (8*eta + 4)*X +  (8*eta - 4)* C*X;
    H = A\B;

    % update Z
    CtC = C'*C;
    A = 2*beta*(2*eta - 1)^2*CtC + mu*I;
    B = 2*beta*(2*eta -1)*(3 - 2*eta)*CtC + Y1 + mu*C;
    Z = A\B;
    Z = Z - diag(diag(Z));
    Z = 0.5*(Z + Z');
    Z = max(Z, 0);
    

    leq1 = C - Z;
    leq2 = C*e - e;

    stopC = max([max(max(abs(leq1))),max(max(abs(leq2)))]);
    
    if iter==1 || mod(iter,50)==0 || stopC<tol
        disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e') ...
            ',leq1=' num2str(max(max(abs(leq1))),'%2.3e')  ...
            ',leq2=' num2str(max(max(abs(leq2))),'%2.3e')]);
    end
    if stopC<tol 
        break;
    else
        Y1 = Y1 + mu*leq1;
        Y2 = Y2 + mu*leq2;
        mu = min(max_mu,mu*rho);
    end

end                                                                                                               

