% BS432
A=[-1  0  0  0;
   0  -12  0  0;
   0  0 -3  0;
   0  0  0 -4];
B=[1  0  0  1;
   0  1  0  1;
   0  0  1  1]';
C=[1  1  0  0;
   0  0  1  1];

K=zeros(3,2);

Kf=[ 2 2;
     2 2;
     0 2];

%K=[-7.2   -7.2/34;
%    14     14/34;
%    0          0];

Q = diag([1 1 1 1]);
V = Q;
R = diag([1 1 1]);

%eig(A-B*K*C)

%a(1)=-7.2;a(2)=14;a(3)=0;a(4)=0.2;

%K=[ a(1)   a(1)*(a(4)-a(3))/(a(1)+a(2));
%    a(2)   a(2)*(a(4)-a(3))/(a(1)+a(2));
%    a(3)*(a(1)+a(2))/(a(4)-a(3))  a(3)]

eig(A-B*K*C)

L =diag([-1 -3 -4 -12]);
%L =diag([-1 -3 -6 -11]);

