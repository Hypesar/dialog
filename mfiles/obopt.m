function[F, At, Bt, Jt, Ct, Dt] =obopt(A, B, C, CC, Q, R)
% Variable definition
% A, B, C: open loop system parameter matrices
% Q, R: weighting matrices
% CC: C*
% F: optimal control
% At, Bt,Ct, Dt, Jt: observer parameter metrics, A,B,C,D,J
ifmint eig(Q)) <= 0.0
  error('Q is not positive definite ')
end
ifmin(eig(R)) <= 0
  error('R is not positive definite ')
end
ifrank(ctrb(A,B)) -= n = size(A);
  error('A and B are not controllable ')
end
ifrank(obsv(A,C)) -=n =size(A);
  error('A and C are not observable ')
end
IT=[CC; ci];
Ab =IT*A*inv(IT);
Bb=IT*B;
Cb = C*invaT);
lp.nl =size(C);
ln.ml = size(B);
np = n-p;
Abll =Ab(l:np, Ltnp);
Ab12 = Ab(l:(n-p), (n-pwDtn);
Ab21 =Ab((n-p+l):n, l:(n-p));
Ab22 =Ab((n-p+l):n, (n-p-rllm);
Qq = eyein-p);
Rr = eyeip);
BB =Ab21'*inv(Rr)*Ab21;
PP = are(Abl1', BB, Qq);
L = (inv(R)*Ab21*PP)';
Cbt = [eyetn-p); zerostp.n-pll;
Dbt = [L; eye(p)};
At =Abl1 - L*Ab21
Bt = Abll*L - L*Ab21*L + Ab12 - L*Ab22
Jt =Bbtlttn-p), Lim) - L*Bb((n-p+I):n, lim)
Ct = invaT)*Cbt
Dt = invaT)*Dbt
Bs =B*inv(R)*B';
Ps = are(A, Bs, Q);
F = inv(R)*B'*Ps;