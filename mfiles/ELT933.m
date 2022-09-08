%
A=[-2.4449   -0.1600         0   -0.1624         0         0   -0.1296         0         0;
    2.5571         0         0   -0.0074         0         0   -0.0059         0         0;
    3.0600         0   -2.0000    0.9843         0         0    0.7854         0         0;
   -0.1624         0         0   -2.5459   -0.1600         0   -0.1111         0         0;
   -0.0074         0         0    2.5525         0         0   -0.0050         0         0;
    0.9843         0         0    3.6720         0   -2.0000    0.6732         0         0;
   -0.1296         0         0   -0.1111         0         0   -2.2892   -0.1600         0;
   -0.0059         0         0   -0.0050         0         0    2.5641         0         0;
    0.7854         0         0    0.6732         0         0    2.1165         0   -2.0000];
%
B=[ 0.9000         0         0;
   -1.0000         0         0;
         0         0         0;
         0    0.9000         0;
         0   -1.0000         0;
         0         0         0;
         0         0    0.9000;
         0         0   -1.0000;
         0         0         0];
%
C=[  0     0     1     0     0     0     0     0     0;
     0     0     0     0     0     1     0     0     0;
     0     0     0     0     0     0     0     0     1];
%
% Lösungen optimale Ausgangsrückführung
K=[ 0.2544    0.0534    0.0402;
    0.0534    0.2864    0.0360;
    0.0402    0.0360    0.2069];
%
%Kf=[ 0 0 2;
%     0 2 0;
%     2 0 0];
%
% K= [ 0.0000    0.0000    0.1041;
%      0.0000    0.3084    0.0000;
%      0.0834    0.0000    0.0000];
%
%K =[ 0.2750         0         0;
%         0    0.3041         0;
%         0         0    0.2230]; 
%
R=diag([ 1 1 1]);
%
Q=diag([ 1 1 1 1 1 1 1 1 1]);
%
V=diag([ 1 1 1 1 1 1 1 1 1]);
%	