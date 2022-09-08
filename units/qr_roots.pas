procedure  qr_root(N:Integer;A : TArrMatReal; Var EwOut,EwIn: TArrVekReal;
                 Var WEnd:Stabil: Boolean; Epsilon : Real); 

(* Eigenwerte aus der Diagonalen selektieren *)
Var
   i,j             : Integer;
   p,q,D,ieps      :  Real;
   b11,b12,b21,b22 : Real;
  // Stabil,Ende     : Boolean;
   
   
Begin 
 
i:=0; j:=0;ieps:=1e-3;
while i< n do begin
      EwOut[j,0]:=A[i,i];   EwOut[j,1]:=0; 
      if i<n then begin
        if abs(A[i+1,i]) > ieps then begin
          b11:=A[i,i];    b12:=A[i,i+1]; 
          b21:=A[i+1,i];  b22:=A[i+1,i+1];          
          p:=-(b11+b22)/2;  q:=b11*b22-b12*b21; D:=p*p - q;          
          if D<0 then begin
             EwOut[j,0]     :=-p; EwOut[j,1]:=sqrt(abs(p*p-q));      //Real+Imag
             EwOut[j+1,0]   :=-p; EwOut[j+1,1]:=-sqrt(abs(p*p-q));   //Real-Imag 
          else 
             EwOut[j,0]     := - p + sqrt(abs(p*p-q));               //Real
             EwOut[j+1,0]   := - p - sqrt(abs(p*p-q));               //Real 
          end; 
          i:=i+1; j:=i;		  
        end;
      end;
	  
	 (* Grenzfälle 0 , N-1 Sonderbehandlung *) 
     if (i=n-1) and (i>0) then begin         
         if abs(A(i,i-1]) > ieps  then begin            
           b11:=A[i-1,i-1];   b12:=A[i-1,i]; 
           b21:=A[i,i-1];     b22:=A[i,i];
           p:=-(b11+b22)/2;     q:=b11*b22-b12*b21; D:=p*p - q;
           j:=n-1; 		  
           if D<0 then begin
              EwOut[j,0]     := -p;   EwOut[j,1]  :=+sqrt(abs(p*p-q));
              EwOut[j+1,0]   := -p;   EwOut[j+1,1]:=-sqrt(abs(p*p-q));            
           else 
              EwOut[j,0]     := - p + sqrt(abs(p*p-q));
              EwOut[j+1,0]   := - p - sqrt(abs(p*p-q));
           end;  	  
         end;
     end;
     i:=i+1; 
     j:=i;
end; // while i<= n do begin
(* norm2 und Stabil berechnen*)
Stabil:=true;seo:=0;sei:=0;
for i:=0to n-1 do begin
    if EwOut[i,1]>0 then begin
       Stabil:=false;   // Instavil
    end;
	seo:=seo + EwOut[i,1]*EwOut[i,1];
    sei:=sei + EwOut[i,1]*EwOut[i,1];	
end;

seo:=sqrt(seo);
sei:=sqrt(sei);
Ende:=true;
abr:=abs(seo-sei);
if abr < Epsilon then  WEnd := false;    

end.