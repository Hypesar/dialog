% Boris Lohmann at - Automatisierungstechnik 39 (1991) 10 Oldenburg Verlag
A=[ -0.0345  -2.23     0      0     2.23;
    47.1    0      -23.55 -23.55  0;
     0      0.304   -0.01   0     0;
     0      0.304    0     -0.01  0;
     0      0        0      0    -8.7];
% 
B=[  0      0        0;
     0      0        0;
     0     -0.608    0;
     0      0       -0.608;
     8.7    0        0];
% 
C=[  0      0.4      0      0     0;
     0      0        1      0     0;
     0      0        0      1     0];
%
K=[ -1.2701    0.1334    0.1552;
     0.0780   -0.7528   -0.0375;
     0.0457   -0.0356   -0.577];
%
V=[ 1 0 0 0 0;
     0 1 0 0 0;
     0 0 1 0 0;
     0 0 0 1 0;
     0 0 0 0 1];
%
Q=[ 1 0 0 0 0;
     0 1 0 0 0;
     0 0 1 0 0;
     0 0 0 1 0;
     0 0 0 0 1];
%
R=[ 1 0 0;
    0 1 0;
    0 0 1];
%