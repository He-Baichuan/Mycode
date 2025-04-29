function [ZZ,TT,SS,alpha,beta,QQ] = schurg(A,B)

    tol = 1e-9;
    [S,T,Q,Z] = qz(A,B);
    S(abs(S) < tol) = 0;
    T(abs(T) < tol) = 0;
    Q(abs(Q) < tol) = 0;
    Z(abs(Z) < tol) = 0;
    [SS,TT,QQ,ZZ] = ordqz(S,T,Q,Z,"udi");
    alpha = diag(SS);
    beta = diag(TT);