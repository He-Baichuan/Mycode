function [C, D, N, L, eigval] = BKsolve(A, B, G, m)

[ZZ,TT,SS,alpha,beta,QQ] =schurg(A,B);
eigval = alpha./beta;
Zp = ZZ';
[a,b] = size(ZZ);

N = inv(Zp(a-m+1:a,a-m+1:a))*Zp(a-m+1:a,1:a-m);
BBN = inv(B(1:a-m,1:a-m)-B(1:a-m,a-m+1:a)*N);
C = BBN*(A(1:a-m,1:a-m)-A(1:a-m,a-m+1:a)*N);
L1 = inv(Zp(a-m+1:a,a-m+1:a))*inv(SS(a-m+1:a,a-m+1:a));
L2 = QQ(a-m+1:a,1:a-m)*G(1:a-m,1)+QQ(a-m+1:a,a-m+1:a)*G(a-m+1:a,1);
L = L1*L2;
D = BBN*(G(1:a-m,1)-A(1:a-m,a-m+1:a)*L);