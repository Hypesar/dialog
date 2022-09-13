% Mehrgrößenreglerentwurf im Zustandsraum
% Prof. Dr.-Ing. U. Konigorski
% Script - Wintersemester 2014/2015
% Simulation 3-Tanksystem Kap. 3.5 S.22
% https://www.arndt-bruenner.de/mathe/scripts/eigenwert2.htm
%T=[ 0,761988198540932; -1,6007361691096187; 1; % zum Eigenwert -6,201472338219242:
%   -0,8345395441411443; 0,2274520456117616; 1; %  zum Eigenwert -2,545095908776476:
%    1,5725513456002111; 1,3732841234978586; 1]; % zum Eigenwert -0,25343175300428267:
%
A=[  1  2  3;
     4  5  6;
     7  8 9];
%    
B=[ 0.5    0;
    0      0;
    0   0.5];
%   
C=[ 1  0   0;
    0  1   0;
    0  0   1];
%   
%K  = zeros(2,3);
%
%K = [
%    4.0000    4.0000   -0.0000 ;
%   -0.0000    4.0000    2.0000 ];
%
Q=diag([ 1 1 1]);
%
V=diag([ 1 1 1]);
%
R=diag([ 1 1]);
%
