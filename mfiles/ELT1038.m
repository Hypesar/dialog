A=[ -0.1088    0.1046         0    0.0523    0.1046    0.1046    0.1046         0    0.1046    0.1046;
   -0.0843   -0.1667         0         0         0         0         0         0         0         0;
   -0.0662         0         0   -0.1250   -0.2500         0         0         0         0         0;
   -0.3313         0    1.2500   -1.8750   -1.2500         0         0         0         0         0;
         0         0         0    0.0625   -0.1250         0         0         0         0         0;
   -0.0021         0         0         0         0   -0.1000         0         0         0         0;
   -0.0015         0         0         0         0         0         0    1.0000         0         0;
         0         0         0         0         0         0   -0.0024   -0.0588         0         0;
   -0.2250         0         0         0         0         0         0         0   -0.2000         0;
   -0.1520         0         0         0         0         0         0         0         0   -0.1800];
%
B=[      0         0         0;
         0         0         0;
    0.1928         0         0;
    0.9628         0         0;
         0         0         0;
    0.0067         0         0;
    0.0048         0         0;
    0.0001         0         0;
         0    0.2000         0;
         0         0    0.1800];
%
C=[ 1.0000         0         0         0         0         0         0         0         0         0;
      0.0350    0.6580         0    0.3290    0.6580    0.6580    0.6580         0   -0.3420   -0.3420;
     -0.0030   -0.4100         0   -0.2050   -0.4100   -0.4100   -0.4100         0    0.5900   -0.4100;
     -0.0320   -0.2480         0   -0.1240   -0.2480   -0.2480   -0.2480         0   -0.2480    0.7520;
         0    1.0000         0         0         0         0         0         0         0         0;
         0         0         0    0.5000    1.0000         0         0         0         0         0;
         0         0         0         0         0    1.0000         0         0         0         0;
         0         0         0         0         0         0    1.0000         0         0         0];
%
% L sung optimale Ausgangsr ckf hrung
% 30 Iterationen 1E-8)
K=[ 0.0179   -0.8424    0.4294    0.4130    1.1096    1.6607    1.1237    1.2724;
   -0.0007   -0.6415    0.5254    0.1161    0.7480    0.7855    0.7677    0.8118;
    0.0571   -0.5723    0.0773    0.4950    0.6375    0.6763    0.6527    0.7052];
%
%K=[  0   -0.2698   -1.3133         0         0    1.5615   -0.0312         0;
%    -0.1411         0         0    0.2786    0.5198         0         0    0.1452;
%          0   -0.5137   -1.1264         0         0    0.5057   -0.0070         0];
%
%K=[ 0   -0.9319   -0.0148         0    1.0000    1.3661    0.8663    1.0000;
%   -0.2825    1.0000    1.0000    0.9099    0.4447         0         0    0.0681;
%         0   -1.3834   -0.4236         0    1.0000    0.9034    0.9102    1.0000];
%
%
% Kf=[2 0 2 0 2 0 2 0;
%    0 2 0 2 0 2 0 2;
%    2 0 2 0 2 0 2 0];
% K=[-0.1917         0   -0.6036         0    0.6496         0   -0.0183         0;
%         0   -0.3398         0   -0.0513         0    0.7001         0    0.2640;
%   -0.1246         0   -0.2489         0    0.4372         0    0.0133         0];
% Kf=[ 0 1 0 1 0 1 0 1;
%     1 0 1 0 1 0 1 0;
%     0 1 0 1 0 1 0 1];
%
% Kf=[0 2 2 0 1 2 2 1;
%    2 1 1 2 2 0 0 2;
%    0 2 2 0 1 2 2 1];
% 
% 24 Iterationen
% K=[  0   -0.2698   -1.3133         0         0    1.5615   -0.0312         0;
%    -0.1411         0         0    0.2786    0.5198         0         0    0.1452;
%          0   -0.5137   -1.1264         0         0    0.5057   -0.0070         0];
%
% Kf=[2 2 2 2 2 2 2 2;
%    2 2 2 2 2 2 2 2;
%    2 2 2 2 2 2 2 2];
%
R=diag([ 1 1 1]);
%
Q=diag([ 1 1 1 1 1 1 1 1 1 1]);
%
V=diag([ 1 1 1 1 1 1 1 1 1 1]);
%