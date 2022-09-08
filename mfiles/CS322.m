% Christian Sonntag
A=[ -4    0    -3;
     1.5 -2     0.1;
     1    0    -1.5];
B=[  1.25 0;
     0    2;
     0.5  0.1];
C=[ 10    0     0;
     0    0.25  0.3];
 
K=zeros(2,2);
Kf=[ 2 0;0 2];
V=eye(3,3);
R=eye(2,2);
Q=V;
eig(A)
L=diag([-4 -2.5 -2]);
%K=lqr(A,B,Q,R)