unit uTools;

// 07.07.2022

{$mode ObjFPC}{$H+}


interface


uses
  Classes, SysUtils, Strings, uDaten;

CONST

MaxMatSize = 2500;
MaxVekSize = 50;
RealSize = 6;
VekDim=16;
MatDim =49;
DUH =49;   {U,H Dimension 7x7 muss Hesse-Dimension einschliessen!}
EPS = 1E-6; {Abruchschranke in VLAK,REGULA }

TYPE
RealMat  = ARRAY[1..MaxMatSize] OF REAL;
RealVek  = ARRAY[1..MaxVekSize] OF REAL;
IntVek   = ARRAY[1..MaxVekSize] of Integer;
Vek      = ARRAY[1..VekDim] OF REAL;
Mat      = ARRAY[1..MatDim] OF REAL;
MUH      = ARRAY[1..DUH] OF REAL;

TarrayReal = Array of Array of Real;

procedure loadfile ( Var FileName: String; Var A_io,B_io,C_io,K_io,Q_io,R_io,V_io; Var n,m,p : Integer);

Procedure TMinv(n:Integer; var a:TArrayReal;var det:Real;var Richtig:Boolean);

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


procedure loadfile ( Var FileName: String; Var A_io,B_io,C_io,K_io,Q_io,R_io,V_io; Var n,m,p : Integer);
Type
   TMatrixReal = Array of Array of Real;
   TVektorReal = Array of Real;

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
    //  if z=8 then begin
     //   z:=z;
     //   end;
      while (TZ[m] = ' ')  do m:=m+1;
      j:=m;
      while (TZ[j]<> ' ') and  (j<lend) do j:=j+1;
      WD:= Copy(TZ,m,j-m);
      zahl:=StrToFloat(WD);
      Ma[s]:=zahl;   // erst Zeile f??llen
      s := s ;
      m:= j;
    end; // while j < lend do begin
      if z=1 then  begin
         Sdim:=s;      // Anzahl Spalten am Ende
      end;
      z:=z+1;
      if l<>0 Then Begin
        z:=z-1;
        Zdim:=z;      // Zeilenende
        Status:=' ';
      end;
      // diag([]) erkennen umwandeln in Matrix
      diag:=Pos('d',TZ);
      if (diag<>0) and (lend=klzu) then begin
        m :=Sdim;Pnt:=1;
        for i:=1 to m do begin
            zahl:=Ma[i];
           for j:=1 to m do begin
              if i=j then begin
                Ma[Pnt]:=zahl;
              end;
              if i>1 then Ma[i]:=0;
              Pnt:=Pnt+1;
           end;
        end;
      end;
      // Umwandlung beendet

  end;  // end pearse


begin
Sl := TStringList.Create;
Sl.LoadFromFile(FileName);
S1Anzahl:= Sl.Count-1;
Status:=' ';

SetLength(Ax,100);
SetLength(Bx,100);
SetLength(Cx,100);
SetLength(Kx,100);
SetLength(Qx,100);
SetLength(Rx,100);
SetLength(Vx,100);

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

 // R??cktransfer
 for i:=0 to N-1 do begin
    for j:=0 to N-1 do begin
        A[i+1+N*j]:= Ax[i+N*j];
        Q[i+1+N*j]:= Qx[i+N*j];
        V[i+1+N*j]:= Vx[i+N*j];
    end;
 end;

 N:=N;

end;  // loadfile





end.
(*--------------------------------------------------------------------------*)
