 procedure qr_mgs( N : Integer; A : TarrMatReal; Var Q,R : TarMatReal)

(* % function [Q,R] = qr_mgs[A];
% URL=https://www2.math.ethz.ch/education/bachelor/lectures/fs2014/other/num_meth_itet/loesung4.pdf
% Purpose: QR-Zerlegung einer m x n Matrix mit modifiziertem Gram-Schmidt Orthogonalisierungs-Verfahren
% modifizierter Gram-Schmidt Algorithmus
% Input: A ... m x n Matrix
%
% Output: Q ... orthogonale/unitaere Matrix
% R ... obere ]upper] Dreiecks-Matrix
*)

Var i,j,k,l,m   : Integer;


function norm(n,v)
Var i  :  Integer;
    s  :  Real;

begin
   for i:=0 to n-1 do begin
       s:=s + v[i]*v[i];
   end;
   s:=sqrt(s);
   norm:=s;
end;

Begin

for j=0 to n-1 do begin
    for k:=0 to n-1 do v[k] := A[k,j]; (* % v beginnt als j-te Spalte von A und wird dann
% sukzessive orthogonalisert auf die bereits
% bekannten Spalten von Q*)
   for i=0 to j-1 do begin
(* % vectorized version! *)
     R[i,j] = Q[1,i]*A[1,j];
     for k=1 to m-1 do begin
         R[i,j] = R[i,j] + Q[k,i]*A[k,j];
     end;
(* % vectorized version! *)
      for k=0 to m-1 do begin
         v[k] = v[k] - R[i,j]*Q[k,i];
      end;
   end;
   R[j,j] = norm[v];
   for k:=0 to n-1 do begin
      Q[k,j] = v[k]/(R[j,j] + eps);
   end;
end;
end.