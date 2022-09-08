procedure Eigenwerte(N:Integer;A : TArrMatReal; Var Ew : TArrVekReal; Var iteration: Integer; Epsilon : Real;  Var Stabil: Boolean;);

iMax=10000;
Epsilon=1e-6;

Var
   Q : TArrMatReal;
   U : TArrMatReal;
   R : TArrMatReal;
   
 i,j,nok : Integer;
   
    Wend : Boolean; 
   
Begin
 
Setlength(Q,n,n);
Setlength(U,n,n);
Setlength(R,n,n);

for i:=0 to n-1 do begin
    for j:=0 to n-1 begin
	   U[i,j]:=0;	
	   Q[i,j]:=0;
	   R[i,j]:=0;	
	end;
	U[i,i]:=1; Q[i,i]:=1;
	Ew[i,0]:=0;Ew[i,1]:=0;
end;

Wend:=true;i:=0;nok:=0;

while Wend & (i<iMax) do begin
 
// Q=eye(n,n); %Set Q as the identity matrix
    for i:=0 to n-1 do begin
        for j:=0 to n-1 do Q[i,j]:=0;	
        Q[i,i]:=1;
    end;
 
  qr_mgs( N,A,Q,R);
  //U=U*Q;  %U*A1*U'
  MaMul(U , Q, U, N,N,N);
   //  A=R*Q;
  MaMul(R , Q, A, N,N,N);
  qr_root(N,A,Ew,Ew,WEnd,Epsilon); 
  i:=i+1;  
end;

end.