unit uTools;


// 16.07.2022

{$mode ObjFPC}{$H+}


interface


uses
  Classes, SysUtils, Strings, uDaten;

Type

TarrayReal = Array of Array of Real;

procedure loadfile ( Var FileName: String; Var A_io,B_io,C_io,K_io,Q_io,R_io,V_io: Mat; Var n,m,p : Integer);

Procedure TMinv(n:Integer; var a:TArrayReal;var det:Real;var Richtig:Boolean);

Procedure MTRSP (VAR A_io : MAT;N,M : Integer );

Procedure Einselemente_Matrix (Var M_o:Mat;Anz:Integer);

Procedure MAFAK (fak : Real;Var M_IO : MAT;N,M :Integer);

Procedure MASUB (VAR A,B,C:MAT;N,M:Integer);

Procedure MASUH (VAR A_i:MUH;Var B_i:Mat;Var C_o:MUH;N,M:Integer);

PROCEDURE MMAMA (A_io,B_io :Mat; Var C :Mat; N,M,P: INTEGER);

Procedure A_BKC (N,M,P: Integer; A_io,B_io,K_io,C_io, Ag_io: MAT);

Procedure QX_Aus (Var Q,C,K,R,QX: MAT;N,M,P: Integer);

Procedure QY_Aus (Var V,B,K,W,QY: MAT;N,M,P: Integer);

Procedure Mat_zu_Vek (VAR QF,ZF,FM:MAT;Var NK,M,P:Integer);

Procedure Vek_zu_Mat  (VAR QF,ZF,FM:MAT;VAR NK,M,P:Integer);



implementation
Procedure TMinv(n:Integer; var a:TArrayReal;var det:Real;var Richtig:Boolean);
{Inversion der Matrix a(i,j),i=j=1..n,Richtig=True}
type TarrayInteger = array of Integer;

var ix,iy,iq : TarrayInteger;
   max,d,t   : Real;


   i,j,k,l,js,iz,m : Integer;
Begin
   //m:= round((n*n+1)/2);  // Anzahl der Unbekannten
   SetLength(ix,n);
   SetLength(iy,n);
   SetLength(iq,n);

   Richtig:=True;det:=1;
   if n>10 then Richtig:=False
   Else begin
     for j:=0 to n-1 do iq[j]:=0;

     for j:=0 to n-1 do begin
       max:=0;
       for i:=0 to n-1 do begin
         if iq[i]<>1 then begin
           for k:=0 to n-1 do
           if(iq[k]<>1)and(max<=abs(a[i,k]))then begin
             js:=k;iz:=i;
             max:=abs(a[i,k]);

           end;
         end;
       end;
       iq[js]:=iq[js]+1;
   {Pivotierung}
       if iz<>js then for i:=0 to n-1 do begin
         t      :=a[iz,i];
         a[iz,i]:=a[js,i];
         a[js,i]:=t;
       end;
       ix[j]:=iz;iy[j]:=js;
   {Kontrolle auf Boesartigkeit}
       d        :=a[js,js];
       det:=det*d;
       if abs(d)<=1E-6 then Richtig:=False
       else begin
         a[js,js] :=1;
         for i:=0 to n-1 do a[js,i]:=a[js,i]/d;
   {Reduktion aller anderen Zeilen}
         for l:=0 to n-1 do if l<>js then begin
           t      :=a[l,js];
           a[l,js]:=0;
           for i:=0 to n-1 do a[l,i]:=a[l,i]-a[js,i]*t;
         end;
       end;
     end;
     if richtig then begin
   {Ruecktausch}
       for j:=0 to n-1 do begin
         i:=n+1-j;
         if ix[i]<>iy[i] then begin
           det:=det*-1;
           iz:=ix[i];js:=iy[i];
           for l:=0 to n-1 do begin
             t      :=a[l,iz];
             a[l,iz]:=a[l,js];
             a[l,js]:=t;
           end;
         end;
       end;
     end;
   end;
 end;


procedure loadfile ( Var FileName: String; Var A_io,B_io,C_io,K_io,Q_io,R_io,V_io: Mat; Var n,m,p : Integer);
Type
   TMatrixReal = Array of Array of Real;
  // TVektorReal = Array of Real;

   TVektorReal    = Array of Real;

Var
  A    : Mat Absolute A_io;
  B    : Mat Absolute B_io;
  C    : Mat Absolute C_io;
  K    : Mat Absolute K_io;
  Q    : Mat Absolute Q_io;
  R    : Mat Absolute R_io;
  V    : Mat Absolute V_io;



 Sl    : TStringList;
 S1Anzahl : Integer;
 TZ    : String;
Status : String;
 i,j   : Integer;
 s,z,x   : Integer;
 Pend  : Integer;
 sx    : Real;

 Ax    : TVektorReal;
 Bx    : TVektorReal;
 Cx    : TVektorReal;
 Kx    : TVektorReal;
 Qx    : TVektorReal;
 Rx    : TVektorReal;
 Vx    : TVektorReal;




  Procedure nuldiag(Var Ma :TVektorReal; dim : Integer);
  var


   i,j   : Integer;

    begin
       for i:=1 to dim do  for j:=1 to dim do  if i<>j then Ma[i]:=0.0;
    end;



  Procedure pearse(var TZ,Status:String;Kz:String;var Ma : TVektorReal; var z,s,Zdim,Sdim : Integer);
  var
    //  Ma    : RealMat Absolute Mat_io;
      zahl     : real;
      WD : String;
     i,j,l,m,lend,klzu, diag    : Integer;
     Pnt : Integer;

  begin
    Status:=Kz;
  // bis [ ]; lesen

     m := Pos('[',TZ);l:=Pos(']',TZ); lend:=Pos(';',TZ);klzu:=Pos(']',TZ);
     if (z=1) and (klzu<>0) then begin
        lend:=klzu; // nur eine Zeile
     end;
  // '.' in ',' umwandeln
     i:=1;
     while (lend>1) and (i<lend) do begin
         if TZ[i] = '.' then TZ[i]:=',';
         i:=i+1;
     end; // (while lend>1) and (i<lend)
     if (l<>0) and (l < lend) then lend:=l;
     m:=m+1; j:=1;
    while j < lend do begin
      if z=8 then begin
       z:=z;
      end;
      while (TZ[m] = ' ')  do m:=m+1;
      j:=m;
      while (TZ[j]<> ' ') and  (j<lend) do j:=j+1;
      WD:= Copy(TZ,m,j-m);
      zahl:=StrToFloat(WD);
      Ma[s]:=zahl;   // erst Zeile füllen
      s := s +1;
      m:= j;
    end; // while j < lend do begin
      if z=1 then  begin
         Sdim:=s;      // Anzahl Spalten am Ende
      end;
      z:=z+1;
      if l<>0 Then Begin
        z:=z-1;
        Zdim:=z;      // Zeilenende
        if Status='R' then begin
          z:=z;
        end;
        Status:=' ';
      end;
      // diag([]) erkennen umwandeln in Matrix
      diag:=Pos('d',TZ);
      if (diag<>0) and (lend=klzu) then begin
        m :=Sdim;Pnt:=0;
        for i:=0 to m-1 do begin
           zahl:=Ma[i];
           for j:=0 to m-1 do begin
              if i=j then begin
                  Ma[Pnt]:=zahl;
                end else begin
              //    Ma[Pnt]:=0;
              end;
              Pnt:=Pnt+1;
           end;
           if (i>0) AND (i<=m) then Ma[i]:=0;    // Rest der Zeile 0 setzen
        end;
        z:=z;
      end;
      // Umwandlung beendet
  end;  // end pearse


begin

  SetLength(Ax,1000);
  SetLength(Bx,1000);
  SetLength(Cx,1000);
  SetLength(Kx,1000);
  SetLength(Qx,1000);
  SetLength(Vx,1000);
  SetLength(Rx,1000);

Sl := TStringList.Create;
Sl.LoadFromFile(FileName);
S1Anzahl:= Sl.Count-1;
Status:=' ';


 for i:=0 to Sl.Count-1 do begin



    TZ := Sl[i];            // ganze Zeile zuweisen

    case TZ[1] of
    '%' : Begin End;       // Kommentarzeile
    'A' : Begin
          z:=1;s:=0;
          pearse(TZ,Status,'A',Ax,z,s,x,N);
          End;
    'B' : Begin
          z:=1;s:=0;
          pearse(TZ,Status,'B',Bx,z,s,x,M);
          End;
    'C' : Begin
          z:=1;s:=0;
          pearse(TZ,Status,'C',Cx,z,s,x,x);
          End;
    'K' : Begin
          z:=1;s:=0;
          pearse(TZ,Status,'K',Kx,z,s,x,x);
          End;
    // Q,R oder V als Diagonaldefinition
    'Q' : Begin
          z:=1;s:=0;
          pearse(TZ,Status,'Q',Qx,z,s,x,N);
     //     nuldiag(Q,N);
          End;
    'R' : Begin
          z:=1;s:=0;
          pearse(TZ,Status,'R',Rx,z,s,x,M);
      //    nuldiag(R,M);
          End;
    'V' : Begin
          z:=1;s:=0;
          pearse(TZ,Status,'V',Vx,z,s,x,N);
      //    nuldiag(V,N);
          End;

     else Begin
        // Folgezeilen zur Matrix
        case Status of
                   ' ' : Begin End;    // Leerzeile
                   'A' : Begin
                         pearse(TZ,Status,'A',Ax,z,s,N,x);
                         End;
                   'B' : Begin
                         pearse(TZ,Status,'B',Bx,z,s,x,x);
                         End;
                   'C' : Begin
                        // s:=1;
                         pearse(TZ,Status,'C',Cx,z,s,P,x);
                         End;
                   'K' : Begin
                        // s:=1;
                         pearse(TZ,Status,'K',Kx,z,s,x,x);
                         End;
                   'Q' : Begin
                         //s:=1;
                         pearse(TZ,Status,'Q',Qx,z,s,x,x);
                         End;
                   'R' : Begin
                         //s:=1;
                         pearse(TZ,Status,'R',Rx,z,s,x,x);
                         End;
                   'V' : Begin
                         //s:=1;
                         pearse(TZ,Status,'V',Vx,z,s,x,x);
                         End;
         end; // case Status of

     end; // else Begin

   end;  // case TZ[1] of

 end; // for i:=0 to Sl.Count do begin

Sl.Free;

// Rücktransfer
for i:=0 to N-1 do begin
   for j:=0 to N-1 do begin

       A[j+1+N*i]:= Ax[j+N*i];
       Q[j+1+N*i]:= Qx[j+N*i];
       V[j+1+N*i]:= Vx[j+N*i];
   end;
end;

for i:=0 to N*N-1 do begin
       A[i+1]:= Ax[i];
       Q[i+1]:= Qx[i];
       V[i+1]:= Vx[i];
end;

for i:=0 to M-1 do begin
   for j:=0 to M-1 do begin
       R[i+1+M*j]:= Rx[i+M*j];
   end;
end;

for i:=0 to N-1 do begin
   for j:=0 to M-1 do begin
       B[j+1+M*i]:= Bx[j+M*i];
   end;
end;

for i:=0 to N-1 do begin
   for j:=0 to P-1 do begin
       C[i+1+N*j]:= Cx[i+N*j];
   end;
end;

for i:=0 to M-1 do begin
   for j:=0 to P-1 do begin
       K[j+1+P*i]:= Kx[j+P*i];
   end;
end;

n:=n;

end;  // loadfile


(*-----------------------------------------------------------*)

{Matrizen - Arithmetik}
 Procedure MTRSP (VAR A_io : MAT;N,M : Integer );
{ Matrix-Transsponierung vom Typ (N,M) }
var i,j : Integer;
    tm  : Mat;
begin
    TM:=A_io;
    for i:=1 to n do for j:=1 to m do
       A_io[j+Pred(i)*M]:=TM[i+Pred(j)*N];
end;

Procedure Einselemente_Matrix (Var M_o:Mat;Anz:Integer);
var j:Integer;
Begin
  For j:=1 to Anz do M_o[j]:=1;
End;

Procedure MAFAK (fak : Real;Var M_IO : MAT;N,M :Integer);
var j:Integer;
Begin
  For j:=1 to N*M do M_IO[j]:=fak*M_IO[j];
End;


Procedure MASUB (VAR A,B,C:MAT;N,M:Integer);
var j:Integer;
Begin
  For j:=1 to N*M do C[j]:=A[j]-B[j];
End;



  Procedure MMAMA (A_io,B_io:Mat; Var C :Mat; N,M,P: INTEGER);


  Var
     A: Mat ABSOLUTE A_io;
     B: Mat ABSOLUTE B_io;
  //   C: Mat ;

     Ax : TarrayReal;
     Bx : TarrayReal;
     Cx : TarrayReal;



     i,j,k : Integer;
     s,x,y     : Real;


     begin

        SetLength(Ax,N,M);
        SetLength(Bx,M,P);
        SetLength(Cx,N,P);

      // Transfer  A,B -> Ax,Bx
      for i:=0 to N-1 do
          for j:=0 to M-1  do
              Ax[i,j]:=A[j+1+(i)*M];

      for i:=0 to M-1 do
          for j:=0 to P-1  do
              Bx[i,j]:=B[j+1+(i)*P];


        for i:=0 to N-1 do begin              // A-Zeilen
          for j:=0 To P-1 do begin            // B-Spalten
                 s:=0;

                 for k:=0 to M-1 do begin
                   x:=Ax[i,k];
                   y:=Bx[k,j];
                   s:=s + x*y;
                 //  s:=s+A[i+(k-1)*N]* B[k+(j-1)*M];  //A-Spalten
                 end;
                 Cx[i,j]:=s;
          end;
          s:=s;
        end;
       n:=n;

      // Rücktransformation

      for i:=0 to N-1 do
           for j:=0 to P-1 do C[j+1+ P*i]:= Cx[i,j];

      Finalize(Ax);
      Finalize(Bx);
      Finalize(Cx);



     end;

 (*
  PROCEDURE MMAMA (A_in,B_in,C_out:Mat; N,M,P: INTEGER);

  VAR
      A: RealMat ABSOLUTE A_in;
      B: RealMat ABSOLUTE B_in;
      C: RealMat ABSOLUTE C_out;

 //  C1: ^Datfeld;
    i,j,k,j1: INTEGER;
    matsize: LongInt;
    s: Extended;
BEGIN

//  getmem(C1,matsize);

  FOR i:=0 TO N-1 DO
    FOR j:=0 TO L-1 DO
      BEGIN
        j1:=Pred(j)*M;
        s:=0.0;
        FOR k:=0 TO M-1 DO s:=s+A[i+k*N]*B[k+j*M];
        C[i+(j)*M]:=s
      END;

END;   { MMAMA }
*)

Procedure MASUH (VAR A_i:MUH;Var B_i:Mat;Var C_o:MUH;N,M:Integer);
var j:Integer;
Begin
  For j:=1 to N*M do C_O[j]:=A_I[j]-B_I[j];
End;

Procedure A_BKC (N,M,P: Integer; A_io,B_io,K_io,C_io, Ag_io: MAT);




Var  i: Integer;
  A    : Mat Absolute A_io;
  B    : Mat Absolute B_io;
  C    : Mat Absolute C_io;
  K    : Mat Absolute K_io;
 // Ag    : Mat Absolute Ag_io;
  Hm,Hz   : Mat;


Begin
  MMAMA(B,K,Hm,n,m,p);    // B*K n,p * p,m => (n,p)
  MMAMA(Hm,C,Hz,n,p,n);  // B*K*C n,p * p,n=> (n,n)
  MASUB (A,Hz,Ag,n,n);   // A-B*K*C
End; { A_BKC }

Procedure QX_Aus (Var Q,C,K,R,QX: MAT;N,M,P: Integer);
Var i,j: Integer;
    test: Boolean;
Begin
  With ARBMAT^ do begin
     { QX := Q + (KC)'RKC }
       MMAMA (K,C,Hm,M,P,N);
       MTRSP (Hm,M,N);
       MMAMA (Hm,R,QX,N,M,M);
       MMAMA (QX,K,Hm,N,M,P);
       MMAMA (Hm,C,QX,N,P,N);
       For i:=1 TO n*n Do QX[i]:=Q[i]+QX[i];
   {    Maus2D  (Qx,n,n,'Qx');}
  End;
End;

Procedure QY_Aus (Var V,B,K,W,QY: MAT;N,M,P: Integer);
Var i,j: Integer;
   Test: Boolean;
Begin
  With ARBMAT^ do begin

 { QY:= V + BKWK'B' statisch}
      MMAMA  (B,K,Hm,N,M,P);
      MTRSP  (Hm,N,P);
      MMAMA  (W,Hm,QY,P,P,N);
      MMAMA  (K,QY,Hm,M,P,N);
      MMAMA  (B,Hm,QY,N,M,N);
      For i:=1 TO n*n Do QY[i]:=V[i]+QY[i];
{      Maus2D  (Qy,n,n,'Qy');}
  END;
End;

Procedure Mat_zu_Vek (VAR QF,ZF,FM:MAT;Var NK,M,P:Integer);
VAR I,J,K :Integer;
BEGIN { Ausfiltern der freien Reglerparameter }
  Nk:=0;
    for j:=1 to P do
      for i:=1 to M do
        if FM[i+Pred(j)*M]=1 then  { nur die Freien }
          begin
            Inc(Nk);
            ZF[Nk]:=QF[i+Pred(j)*M];
          end;
END;

procedure Vek_zu_Mat  (VAR QF,ZF,FM:MAT;VAR NK,M,P:Integer);
var    k,i,j,ij:   integer;
begin { Umspeichern Ky_Vek ---> Ky_Mat anhand von Ky_LOG }
  Nk:=0;
    for j:=1 to P do
      for i:=1 to M do
      begin
        ij:=i+Pred(j)*M;
        if FM[ij]=1 then
          begin
            INC(Nk);
            ZF[ij]:=QF[Nk];
          end
      end
end;



end.
(*--------------------------------------------------------------------------*)

