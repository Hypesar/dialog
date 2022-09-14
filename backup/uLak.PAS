Unit uLak;
// 04.07.2022
Interface

uses
(*  Matar,Matgl; *)
 uDaten, uTools;


{--------------------------------------------------------------------}
procedure lyap( N : Integer; Var A_io,Q_io,P_io);

PROCEDURE Ljapunow (N: INTEGER; VAR A_io,Q_io,P_out,U_io;
                                LjpForm2:BOOLEAN;VAR skip_io: BOOLEAN);

PROCEDURE Liapunow (N: INTEGER; VAR A_io,Q_io,P_out,U_io;
                                LjpForm2:BOOLEAN;VAR skip_io: BOOLEAN);

Procedure Eigenwerte (N:Integer;Var Ag_io;Var Ew_io:Vek;Var TM_io; Var Stabil_io:Boolean);

//Procedure Eigenwerte (N:Integer;Var Ag:Mat;Var Ew:Vek;Var Stabil:Boolean);

procedure BALANCE (n: INTEGER; VAR A_io;
                               VAR low,upp: INTEGER; VAR d_out);


procedure BALBAK (n,low,upp,m: INTEGER; VAR d_in, Z_io);



function CABS (xr, xi: REAL): REAL;

PROCEDURE CDIV (xr,xi, yr,yi: REAL; VAR zr,zi: REAL);

function SIGN (x: REAL): INTEGER;

function SPUR (VAR A_in; N: INTEGER): REAL;


PROCEDURE TRANSP (n: INTEGER; VAR A_io);

PROCEDURE ORTHES (n,k,l: INTEGER; VAR A_io, d_out);


PROCEDURE ORTRANS (n,low,upp: INTEGER; VAR H_in, V_out, d_io);



PROCEDURE HQR2S (n,low,upp: INTEGER;
                 VAR H_io, Vecs_io, wr_out, wi_out, cnt_out;
                                                Schur: BOOLEAN);


PROCEDURE EINORM (n: INTEGER; VAR Vec_in, Vec_out, wi_in);


PROCEDURE EIGRF (n: INTEGER; VAR A_io,V_out,wr_out,wi_out);


PROCEDURE LIAPS (n: INTEGER; VAR A_in,Q_io,P_out);


PROCEDURE UTCU (n: INTEGER; VAR U_in,C_io; rueck: BOOLEAN);

PROCEDURE ROCOEX (n: INTEGER; VAR A_io);


PROCEDURE Lia_Probe (n: INTEGER; VAR A_in,Q_io,P_in; LjpForm2: BOOLEAN);

{====================================================================}

Implementation

CONST
      MaxMatSize = 2500;
      MaxVekSize = 50;
      RealSize   = 6;       (* LÑnge der Realzahldarstellung *)
      basis      = 2.0;     {Basis der REAL-Zahl-Darstellung}
      eta        = 1.0E-38; {kleinste darstellbare REAL-Zahl}
      macheps    = 2.0E-12; {kleinste pos. Zahl fuer 1+macheps > 1.0}
      Tol        = 2.0E-26; {eta/macheps}
      eps8       = 1.5E-11; { eps8 = 8*macheps  }


TYPE
      RealMat  = ARRAY[1..MaxMatSize] OF REAL;
      RealVek  = ARRAY[1..MaxVekSize] OF REAL;
      IntVek   = ARRAY[1..MaxVekSize] of Integer;


Const
   eps=1.0E-15;    (* fuer MINV1 *)
   eps1=1.0E-7;    (* fuer EXPAT *)
   epsx=1.0E-4;    (* fÅr  MAEIG *)
   maxitera=31;    (* fuer EXPAT *)
   Feldmax=2500;   (* max. 50 * 50 Matrizen *)


Type
   Datfeld = Array[1..Feldmax] of real;

procedure lyap(N : Integer; Var A_io,Q_io,P_io);
// Ljapunovgleichung  A*P + P*A' = Q
Type
   // Mat     =  array[1..100] of real;

    TArrayReal    = array of array of Real;
    TArrayInteger = array of array of Integer;
Var
 A    : RealMat absolute A_io;
 Q    : RealMat absolute Q_io;
 P    : RealMat absolute P_io;

// Zwischenspeicher

  u,v,l,k,m,i,j : Integer;
  DetA         : Real;
  Richtig      : Boolean;
//  MaInvSel     : Integer;
 //W,Ux,Vx: ^RealMat;
  Pnt     : TArrayInteger;
  Wx,Wt,Ux,Vx,At,Ax,Mp   : TArrayReal;
  Qp, Xp , Px,Mw,R,Rt         : TArrayReal;


(* Matrix - Invertierung *)


procedure MaInv(N: Integer; Ma : TArrayReal; Var Det : Real; Var Richtig : Boolean);
(* Matrix - Invertierung *)
Type

   TVektorInteger = Array of Integer;
Var
   i,j,k,l,m : Integer;

   pv,Dmax   : Real;
   iD,iM   : TVektorInteger;

 begin
   SetLength(iD,N);
   SetLength(iM,N);
   for i:=0 to N-1 do begin
       iD[i]:=1;  // Abarbeitung Maximum
       iM[i]:=0;
   end;

   Det:=1; k:=0; l:=0;
   Richtig:=True;

 //  for k:=0 to N-1 do
   while (l<N) AND (Richtig=True)do begin
      // Diagonale Maximum suche
      Dmax:=0;
      for i:=0 to N-1 do begin
          if (Dmax < abs(Ma[i,i])) And (id[i]<>0) then begin
             Dmax:=Ma[i,i];
             k:=i;
          end;
      end;
      id[k]:=0;       // aktuelle Diagonalposition austragen
      im[l]:=k;       // Abarbeitungsreihenfolge merken
      for m:=0 to l do begin
          if (im[m]<> k) And (iD[m]<>k) then im[l]:=k;
      end;

      det:=det*Ma[k,k];
      Ma[k,k]:=-Ma[k,k];
      if Ma[k,k]<>0 Then begin
        pv := 1/ Ma[k,k];
         for i := 0 to N-1 do begin
           Ma[k,i]:=Ma[k,i]*pv;
         end;
         for i:=0 to N-1 do begin
           if i<>k Then begin
              for j:=0 to N-1 do Begin
                 if j<>k then begin
                    Ma[i,j]:= Ma[i,j] + Ma[k,j]*Ma[i,k];
                 end;
              end;
           end;
           Ma[i,k] := -Ma[i,k]*pv;
         end;
        // k:=k + 1;
         l:=l + 1;
      end //  if Ma[k,k]<>0 Then begin

      else Richtig:=False;

     end;
     det:=det;
end;

procedure minv(N: Integer; Ma : TArrayReal; Var Det : Real);
Var
   i,j,k : Integer;
   pv   : Real;

 begin
   Det:=1;
   for k:=0 to N-1 do begin
      det:=det*Ma[k,k];
      Ma[k,k]:=-Ma[k,k];
         pv := 1/ Ma[k,k];
      for i := 0 to N-1 do begin
         Ma[k,i]:=Ma[k,i]*pv;
      end;
      for i:=0 to N-1 do begin
        if i<>k Then begin
           for j:=0 to N-1 do Begin
              if j<>k then begin
                 Ma[i,j]:= Ma[i,j] + Ma[k,j]*Ma[i,k];
              end;
           end;
        end;
        Ma[i,k] := -Ma[i,k]*pv;
     end;
   end;
   det:=det;
end;



(*-------------------------------------*)

(* Matrix Transponierung *)
Procedure MTRSP (VAR Ma,Mb : TArrayReal;N,M : Integer );
{ Matrix-Transsponierung vom Typ (N,M) }
var i,j : Integer;

begin
    for i:=0 to n-1do for j:=0 to m-1 do Mb[j,i]:=Ma[i,j];
end;


(*-------------------------------------*)

(* Matrizen multiplikation *)

Procedure MaMul(Var Ma,Mb,Mc : TArrayReal; N,M,P : Integer);
Var
   i,j,k,l : Integer;
   s     : Real;
   begin
      for i:=0 to N-1 do begin              // A-Zeilen
        for j:=0 To P-1 do begin            // B-Spalten
               s:=0;  //s:=Ma[i,0]*Mb[0,j];
               for k:=0 to M-1 do s:=s+Ma[i,k]*Mb[k,j];  //A-Spalten
               Mc[i,j]:=s;
        end;
      end;
   end;

Procedure MaSub(Var Ma,Mb,Mc : TArrayReal; N,M : Integer);
Var
   i,j: Integer;

   begin
       for i:=0 to N-1 do begin              // A-Zeilen
        for j:=0 To M-1 do begin            // B-Spalten
             Mc[i,j]:=Ma[i,j] - Mb[i,j];
        end;
      end;
end;

procedure ljpQuality(Var N,M: Integer; Qp,Wx,Xp : TArrayReal; Var QDif: Real);

Type
    TVektorReal = Array of Real;

Var R,Rt,Qd : TArrayReal;

Begin
   SetLength(R,N*N,1);  // Spaltenvektor
   SetLength(Rt,1,N*N); // Transponierte
   SetLength(Qd,1,1); // Skalar
 // R = Wx*x -Q(:);r'*r      % Berechnungsgenauigkeit
   MaMul(Wx,Xp,R,N*N,M,1); // (N*N,m)*(m,1) => (m,1);
   MaSub(R,Qp,R,N*N,1);
   MTRSP(R,Rt,N*N,1);
   MaMul(Rt,R,Qd,1,N*N,1);
   Qdif:=Qd[0,0];

End;

(*-------------------------------------*)

(*-------------------------------------*)


Begin
(* dynamischen Speicher definieren & abrufen*)
m:=round(n*(n+1)/2);  // Anzahl der Unbekannten  in Px
SetLength(Ax,N,N);
SetLength(At,N,N);
SetLength(Px,N,N);
SetLength(Qp,N*N,1);
SetLength(R,N*N,1);
SetLength(Rt,1,N*N);
SetLength(Xp,M,1);
SetLength(Ux,N*N,M);
SetLength(Vx,N*N,M);
SetLength(Wx,N*N,M);
SetLength(Wt,M,N*N);
SetLength(Mw,M,N*N);
SetLength(Mp,M,M);

SetLength(Pnt,N,N);  // Pointer

(* Konv. Vektor 1D to 2D Matrizen *)

for i:=0 to N-1 do
    for j:=0 to N-1 do begin
        Ax[j,i]:= A[i+1+N*j];
        At[i,j]:= A[i+1+N*j];  // Transponierte
    end;
// Qp = Qx(:) Vektor
for i:=0 to N*N-1 do Qp[i,0]:= -Q[i+1];

(* Index Pnt Pointer init *)
k:=0;
for i:=0 to N-1 do begin
    for j:=i to N-1 do begin
        Pnt[i,j]:=k;
        if i<>j then begin
           Pnt[j,i] := Pnt[i,j];
        end;
         k:=k+1;
    end;
end;

// Gleichungssystem aufbauen

k:=0;
for u:=0 to n-1 do begin
    for v:=0 to  n-1 do begin
        for i:=0 to n-1 do begin
           Ux[k,Pnt[i,v]] := Ax[u,i];    // U = A*P
           Vx[k,Pnt[i,u]] := At[i,v];    // V = P*A'
        end;
        k:=k+1;
    end;
end;

//  Wx := Ux + Vx    (n,m)
for i:=0 to M-1 do
    for j:=0 to N*N-1 do
       Wx[j,i] := Ux[j,i] + Vx[j,i];

// Wt = Wx'  (m,n)  transponiert
for i:=0 to M-1 do
    for j:=0 to N*N-1 do
        Wt[i,j] := Wx[j,i];

// Xp = -inv(W'*W)*W'* Q(:);   (m,1)

// Mp = Wt*Wx  (m,n)*(n,m) =>  (m,m)
MaMul(Wt,Wx,Mp,M,N*N,M);

// Mp Inverse berechnen

case MaInvSel of

1: begin
      minv(M, Mp, DetA);
   end;
2: begin
     MAinv(M,Mp,DetA,Richtig);
   end;
3: begin
     TMinv(M,Mp,DetA,Richtig);
   end

else

end;



//Mw = inv(Mp)*Wt   (m,m)*(m,n) => (m,n)
MaMul(Mp,Wt,Mw,M,M,N*N);
// Xp = Mw*Qp       (m,n) * (n,1) => (m,1)
MaMul(Mw,Qp,Xp,M,N*N,1);

ljpQuality(N,M, Qp, Wx, Xp, QDif);   // Berechnungsg¸te Mat Inv.

// R¸cksortierung von x in P via Px
for i:=0 to n-1 do begin
    for j:=0 to n-1 do begin
        Px[i,j]:=Xp[Pnt[i,j],0];
        if i<>j Then begin
           Px[j,i]:=Px[i,j];
        end;
    end;
end;

// R¸cktransfer von Px zu P

for i:=0 to N-1 do
    for j:=0 to N-1 do P[i+1+ N*j]:= Px[i,j];

// Heap Freigaben

Finalize(Ax);
Finalize(At);
Finalize(Qx);
Finalize(Ux);
Finalize(Vx);
Finalize(Wx);
Finalize(Px);
Finalize(Qp);
Finalize(Mp);
Finalize(Pnt);
Finalize(Xp);

end;




PROCEDURE Liapunow (N: INTEGER; VAR A_io,Q_io,P_out,U_io;
                                LjpForm2:BOOLEAN;VAR skip_io: BOOLEAN);

LABEL fail;
(*skip_Io=TRUE: A hat bereits Schurform und Schurvektoren U liegen vor
  skip_iO=False:A ist instabil
    LjpForm2=True : AP + PA' = -Q
    LjpForm2=FALSE: A'P + PA = -Q *)
VAR
  A: RealMat Absolute A_io;
  Q: RealMat Absolute Q_io;
  P: RealMat Absolute P_out;
  U: RealMat Absolute U_io;
Skip: BOOLEAN Absolute Skip_IO;
  wr,wi,cnt: ^RealVek;
  j,low,upp: INTEGER;
  vecsize: LongInt;
  Tm,D:Real;
BEGIN  (* Ljapunow *)
 (* if severe_error then Exit; *)
  vecsize:=n*RealSize;
(*  if maxavail < 3*vecsize then
    begin
      GL_Fehler(12,severe_error); Exit
    end; *)
  getmem(wr,vecsize); getmem(wi,vecsize);
  getmem(cnt,vecsize);
  IF NOT skip THEN  {zunaechst Schur-Trans.}
    BEGIN
      low:=1; upp:=N;
      ORTHES(N,low,upp,A,wr^);
      ORTRANS(N,low,upp,A,U,wr^);
      HQR2S(N,low,upp,A,U,wr^,wi^,cnt^,TRUE);
      Skip:=True;{Stabilitaetstest}
      Tm:=0;D:=0;
      For j:=low To upp Do begin
      (*    Tm:=Tm+1/Abs(wr^[j]);D:=D+1/Sqrt(SQR(wr^[j])+SQR(wi^[j]));  *)
          If wr^[j]>=0 Then skip:=False;
     (*    gotoxy(4,j+1);write ('Re=',wr^[j]:8:4,' Im=',wi^[j]:8:4);  *)
      End;
     (*    gotoxy(3,1);write (' Tm=',Tm:8:4,' D =',Tm/D:8:4);  *)

      if Not Skip then Goto fail;
   (*   if severe_error then Goto fail;*)
    END; { if }
  IF LjpForm2 THEN { Aschur:= T A'schur T }
  BEGIN TRANSP(N,A); ROCOEX(N,A) END;
  UTCU(N,U,Q,FALSE); IF LjpForm2 THEN ROCOEX(N,Q);
  LIAPS(N,A,Q,P); IF LjpForm2 THEN ROCOEX(N,P);
  UTCU(N,U,P,TRUE);
fail:
  freemem(cnt,vecsize);
  freemem(wi,vecsize);
  freemem(wr,vecsize)
END;  (* Liapunow *)



procedure BALANCE (n: INTEGER; VAR A_io;
                               VAR low,upp: INTEGER; VAR d_out);

LABEL  L1,L2,iteration,L3,L4;

VAR
  a: RealMat Absolute A_io;
  d: RealVek Absolute d_out;
  i,j,k,l,ij,ij1:   INTEGER;
  b2,c,f,g,r,s:     Extended;
  noconv:           BOOLEAN;

  procedure exc (m: INTEGER);
  var
    i: Integer;
  BEGIN
    d[m]:=j;
    IF j <> m THEN
      BEGIN
        ij:=n*Pred(j); ij1:=n*Pred(m);
        FOR i:=1 TO k DO
          BEGIN
            f:=a[i+ij];
            a[i+ij]:=a[i+ij1];
            a[i+ij1]:=f;
          END;
        FOR i:=l TO n DO
          BEGIN
            ij:=n*Pred(i); f:=a[j+ij]; a[j+ij]:=a[m+ij]; a[m+ij]:=f
          END
      END
  END;  { exc }

BEGIN  (* BALANCE *)
  b2:=basis*basis; l:=1; k:=n;

  { search for rows isolating an eigenvalue and push them down }

  j:=K;
  while j>=1 do begin
    r:=0.0;
    for i:=1 to j-1 do
      r:= r + abs(a[j+n*(i-1)]);
    for i:=j+1 to k do
      r:= r + abs(a[j+n*(i-1)]);
    if r = 0 then begin
       exc(k); k:=k-1; j:=k;
       end else
       j:=j-1;
    end;
   r:=r;

  { search for columns isolating an eigenvalue and push them left }
  j:=1;
  while j<=k do begin
  c:=0.0;  ij:=n*Pred(j);
    for i:=1 to j-1 do
      c:= c + abs(a[i+ij]);
    for i:=j+1 to k do
      c:= c + abs(a[i+ij]);
    if (c = 0)  then begin
       exc(l); l:=Succ(l); j:=l;
       end else
       j:=j+1;
    end;
  c:=c;

   { search for rows isolating an eigenvalue and push them down }
(*
  L1:
  FOR j:=k DOWNTO 1 DO
    BEGIN
      r:=0.0;
      FOR i:=1 TO Pred(j) DO begin
        r:=r+ABS(a[j+n*Pred(i)]);
      end;
      FOR i:=Succ(j) TO k DO begin
        r:=r+ABS(a[j+n*Pred(i)]);
      end;
      IF r=0.0 THEN
        BEGIN exc(k); k:=Pred(k); GOTO L1; END
    END;

*)
  { search for columns isolating an eigenvalue and push them left }
(*
L2:
  FOR j:=l TO k DO
    BEGIN
      c:=0.0; ij:=n*Pred(j);
      FOR i:=l TO Pred(j) DO c:=c+ABS(a[i+ij]);
      FOR i:=Succ(j) TO k DO c:=c+ABS(a[i+ij]);
      IF c=0.0 THEN
        BEGIN exc(l); l:=Succ(l); GOTO L2; END
    END;
*)
  { now balance the submatrix in rows l through k }
  low:=l; upp:=k;
  FOR i:=l TO k DO d[i]:=1;
iteration:
  noconv:=FALSE;
  FOR i:=l TO k DO
    BEGIN
      c:=0.0; r:=0.0;
      FOR j:=l TO Pred(i) DO
        BEGIN c:=c+ABS(a[j+n*Pred(i)]); r:=r+ABS(a[i+n*Pred(j)]) END;
      FOR j:=Succ(i) TO k DO
        BEGIN c:=c+ABS(a[j+n*Pred(i)]); r:=r+ABS(a[i+n*Pred(j)]) END;
      g:=r/basis; f:=1.0; s:=c+r;


      while c < g do BEGIN
       f:=f*basis; c:=c*b2;
      end;

      while c >= g do BEGIN
       f:=f/basis; c:=c/b2;
      end;

      g:=r*basis;
 (*
      L3:
            IF c < g THEN BEGIN f:=f*basis; c:=c*b2; GOTO L3; END;
L4:
      IF c >= g THEN BEGIN f:=f/basis; c:=c/b2; GOTO L4; END;

*)

      IF (c+r)/f < 0.95*s THEN
        BEGIN
          g:=1/f; d[i]:=d[i]*f; noconv:=TRUE;
          FOR j:=l TO n DO BEGIN ij:=i+n*Pred(j); a[ij]:=a[ij]*g END;
          ij:=n*Pred(i); FOR j:=1 TO k DO a[j+ij]:=a[j+ij]*f
        END
    END;
  IF noconv THEN GOTO iteration;
END; { BALANCE }


{--------------------------------------------------------------------}


procedure BALBAK (n,low,upp,m: INTEGER; VAR d_in, Z_io);


VAR
   z: RealMat Absolute Z_io;
   d: RealVek Absolute d_in;
   i,j,k,ij:INTEGER;
   s: Extended;

BEGIN
  FOR i:=low TO upp DO
    BEGIN
      s:=d[i];
      { left hand eigenvectors are back transformed,          }
      { if the foregoing statement is replaced by s:=1/d[i]   }
      FOR j:=1 TO m DO
        BEGIN ij:=i+n*Pred(j); z[ij]:=z[ij]*s END
    END {i};
  FOR i:=Pred(low) DOWNTO 1 DO
    BEGIN
      k:=Round(d[i]);
      IF k<>i THEN
        FOR j:=1 TO m DO
          BEGIN
            ij:=i+n*Pred(j);
            s:=z[ij]; z[ij]:=z[k+n*Pred(j)]; z[k+n*Pred(j)]:=s
          END {j}
    END {i: erster Teil};
  FOR i:=Succ(upp) TO n DO
    BEGIN
      k:=Round(d[i]);
      IF k<>i THEN
        FOR j:=1 TO m DO
          BEGIN
            ij:=i+n*Pred(j);
            s:=z[ij]; z[ij]:=z[k+n*Pred(j)]; z[k+n*Pred(j)]:=s
          END {j}
    END {i: zweiter Teil}
END; {balbak}

{--------------------------------------------------------------------}



{--------------------------------------------------------------------}


function SPUR (VAR A_in; N: INTEGER): REAL ;


VAR
  A: Datfeld Absolute A_in;
  R: Extended;
  i: INTEGER;
BEGIN
  R:=0.0;
  FOR i:=1 TO N DO begin

    R:=R+A[i+Pred(i)*N];
    end;

  SPUR:=R;
END;

{--------------------------------------------------------------------}


PROCEDURE TRANSP (n: INTEGER; VAR A_io);

VAR
   A: Datfeld Absolute A_io;
   i,j: INTEGER;
   s: REAL;
BEGIN
  FOR i:=1 TO n DO
    FOR j:=Succ(i) TO n DO
      BEGIN
        s:=A[i+n*Pred(j)];
        A[i+n*Pred(j)]:=A[j+n*Pred(i)];
        A[j+n*Pred(i)]:=s
      END
END;

{--------------------------------------------------------------------}


function CABS (xr, xi: REAL): REAL;

VAR
  h,cxr,cxi: Extended;
BEGIN
  cxr:=ABS(xr); cxi:=ABS(xi);
  IF cxi > cxr THEN BEGIN h:=cxr; cxr:=cxi; cxi:=h; END;
  IF cxi=0 THEN CABS:=cxr
           ELSE CABS:=cxr*Sqrt(1.0+Sqr(cxi/cxr));
END;  {CABS}

{--------------------------------------------------------------------}


PROCEDURE CDIV (xr,xi, yr,yi: REAL; VAR zr,zi: REAL);

VAR
  h,cyr,cyi: Extended;
BEGIN

  if yr=0.0 then if yi=0.0 then
    begin
    (*  MA_Fehler(1,severe_error); Exit  *)
    end;
  IF ABS(yr) > ABS(yi) THEN
    BEGIN
      h:=yi/yr; cyr:=h*yi+yr;
      zr:=(xr+h*xi)/cyr; zi:=(xi-h*xr)/cyr;
    END
                       ELSE
    BEGIN
      h:=yr/yi; cyi:=h*yr+yi;
      zr:=(h*xr+xi)/cyi; zi:=(h*xi-xr)/cyi;
    END;
END;  {CDIV}

{--------------------------------------------------------------------}


function SIGN (x: REAL): INTEGER;

BEGIN  { Signumfunktion }
  IF x > 0.0 THEN Sign:=1
             ELSE
    IF ABS(x) = 0.0 THEN Sign:=0
                    ELSE Sign:=-1
END; { SIGN }

{--------------------------------------------------------------------}


PROCEDURE ORTHES (n,k,l: INTEGER; VAR A_io, d_out);



VAR
    a: RealMat Absolute A_io;
    d: RealVek Absolute d_out;
    i,j,m,la,ij: INTEGER;
    f,g,h: Extended;
    x   : Real;
    Label
      skip;

BEGIN
  la:=Pred(l);
  FOR m:=Succ(k) TO la DO
    BEGIN
      h:=0.0; ij:=n*(m-2);
      FOR i:=l DOWNTO m DO
        BEGIN
          f:=a[i+ij];
          d[i]:=f;
          h:=h+sqr(f);
        END;
      IF h <= Tol THEN BEGIN g:=0.0; GOTO skip; END;
      IF f>=0 THEN g:=-sqrt(h) ELSE g:=sqrt(h);
      h:=h-f*g; d[m]:=f-g;

      {Berechnung von (I-[u*u']/h)*A }
      FOR j:=m TO n DO
        BEGIN
          f:=0.0;
          FOR i:=l DOWNTO m DO begin
            f:=f+d[i]*a[i+n*Pred(j)];
          end;
          f:=f/h;
          FOR i:=m TO l DO
            BEGIN ij:=i+n*Pred(j); a[ij]:=a[ij]-f*d[i] END
        END;

      {Berechnung von (I-[u*u']/h)*A*(I-(u*u')/h)  }
      FOR i:=1 TO l DO
        BEGIN
          f:=0.0;
          FOR j:=l DOWNTO m DO begin
            x:=a[i+n*Pred(j)];
            f:=f+d[j]*a[i+n*Pred(j)];
          end;
          f:=f/h;
          FOR j:=m TO l DO
            BEGIN ij:=i+n*Pred(j); a[ij]:=a[ij]-f*d[j] END
        END;
skip:
        a[m+n*(m-2)]:=g;
  END
 END;  { ORTHES }

{--------------------------------------------------------------------}


PROCEDURE ORTRANS (n,low,upp: INTEGER; VAR H_in, V_out, d_io);


LABEL   skip;
VAR
   h: RealMat Absolute H_in;
   v: RealMat Absolute V_out;
   d: RealVek Absolute d_io;
   i,j,k,m,ij: INTEGER;
   x,y: Extended;

BEGIN
  FOR i:=1 TO n DO
    BEGIN
      FOR j:=1 TO n DO v[i+n*Pred(j)]:=0.0;
      v[i+n*Pred(i)]:=1.0
    END{i};
  FOR k:=upp-2 DOWNTO low DO
    BEGIN
      m:= Succ(k); y:=h[m+n*Pred(k)];
      IF ABS(y)=0.0 THEN GOTO skip;
      y:=y*d[m];
      FOR i:=k+2 TO upp DO d[i]:=h[i+n*Pred(k)];
      FOR j:=m TO upp DO
        BEGIN
          x:=0.0;
          FOR i:=m TO upp DO
            x:=x+d[i]*v[i+n*Pred(j)];
          x:=x/y;
          FOR i:=m TO upp DO
            BEGIN
              ij:=i+n*Pred(j); v[ij]:=v[ij]+x*d[i]
            END {i}
        END {j};
skip:
    END {k}
END; { ortrans }


{--------------------------------------------------------------------}


PROCEDURE HQR2S (n,low,upp: INTEGER;
                 VAR H_io, Vecs_io, wr_out, wi_out, cnt_out;
                                                Schur: BOOLEAN);


LABEL   nextw, nextit, cont1, cont2, cont3, onew, twow, fin, fail;

VAR
   h:    RealMat Absolute H_io;
   vecs: RealMat Absolute Vecs_io;
   wr:   RealVek Absolute wr_out;
   wi:   RealVek Absolute wi_out;
   cnt:  IntVek  Absolute cnt_out;
   i,j,k,l,m,na,its,en,ij: INTEGER;
   p,q,r,s,t,w,x,y,z,ra,sa,vr,vi,norm: Extended;
   notlast: BOOLEAN;

BEGIN
 (* if severe_error then exit; *)
  { isolated roots }
  FOR i:=1 TO Pred(low) DO
    BEGIN wr[i]:=h[i+n*Pred(i)]; wi[i]:=0.0; cnt[i]:=0 END;
  FOR i:=Succ(upp) TO n DO
    BEGIN wr[i]:=h[i+n*Pred(i)]; wi[i]:=0.0; cnt[i]:=0 END;
  en:=upp; t:=0.0;

nextw:
  IF en < low THEN GOTO fin;
  its:=0; na:=Pred(en);
  { look for single small sub-diagonal element }

nextit:
  FOR l:=en DOWNTO Succ(low) DO
    BEGIN
      ij:=l+n*(l-2);
      IF ABS(h[ij]) <= macheps*(ABS(h[Pred(ij)])+ABS(h[ij+n])) THEN
        GOTO cont1;
    END;
  l:=low;

cont1:
  x:=h[en+n*Pred(en)]; IF l = en THEN GOTO onew;
  y:=h[na+n*Pred(na)]; w:=h[en+n*Pred(na)]*h[na+n*Pred(en)];
  IF l = na THEN GOTO twow;
  IF its = 30 THEN
    BEGIN
    (*  GL_Fehler(3,severe_error); cnt[en]:=31; GOTO fail  *)
    END;
  IF (its=10) OR (its=20) THEN
    BEGIN  { form exceptional shift }
      t:=t+x;
      FOR i:=low TO en DO
        BEGIN ij:=i+n*Pred(i); h[ij]:=h[ij]-x END;
      s:=ABS(h[en+n*Pred(na)])+ABS(h[na+n*(en-3)]);
      y:=0.75*s; x:=y; w:=-0.4375*sqr(s)
    END;
  Inc(its);
  { look for two consecutive small sub-diagonal elements }
  FOR m:=en-2 DOWNTO l DO
    BEGIN
      ij:=m+n*Pred(m);
      z:=h[ij]; r:=x-z; s:=y-z;
      p:=(r*s-w)/h[Succ(ij)]+h[ij+n];
      q:=h[Succ(ij+n)]-z-r-s; r:=h[ij+n+2];
      s:=ABS(p)+ABS(q)+ABS(r); p:=p/s; q:=q/s; r:=r/s;
      IF m = l THEN GOTO cont2;
      IF ABS(h[ij-n])*(ABS(q)+ABS(r)) <=
         macheps*ABS(p)*(ABS(h[ij-Succ(n)])+ABS(z)+ABS(h[ij+Succ(n)]))
        THEN GOTO cont2;
    END; { m }

cont2:
  FOR i:=m+2 TO en DO h[i+n*(i-3)]:=0.0;
  FOR i:=m+3 TO en DO h[i+n*(i-4)]:=0.0;
  { double QR step involving rows l to en and columns m to en of the
    complete array }
  FOR k:=m TO na DO
    BEGIN
      notlast:=k <> na; ij:=k+n*Pred(k);
      IF k <> m THEN
        BEGIN
          p:=h[ij-n]; q:=h[ij-Pred(n)];
          IF notlast THEN r:=h[ij-n+2] ELSE r:=0.0;
          x:=ABS(p)+ABS(q)+ABS(r);
          IF x=0.0 THEN GOTO cont3;
          p:=p/x; q:=q/x; r:=r/x
        END;
      s:=Sqrt(sqr(p)+sqr(q)+sqr(r));
      IF p < 0.0 THEN s:=-s;
      IF k <> m THEN h[ij-n]:=-s*x
                ELSE IF l <> m THEN h[ij-n]:=-h[ij-n];
      p:=p+s; x:=p/s; y:=q/s; z:=r/s;
      q:=q/p; r:=r/p;
      { row modification }
      FOR j:=k TO n DO
        BEGIN
          ij:=k+n*Pred(j);
          p:=h[ij]+q*h[Succ(ij)];
          IF notlast THEN
            BEGIN
              p:=p+r*h[ij+2];
              h[ij+2]:=h[ij+2]-p*z
            END;
          h[Succ(ij)]:=h[Succ(ij)]-p*y;
          h[ij]:=h[ij]-p*x
        END;{j}
      IF k+3 < en THEN j:=k+3 ELSE j:=en;
      { column modification }
      FOR i:=1 TO j DO
        BEGIN
          ij:=i+n*k;
          p:=x*h[ij-n]+y*h[ij];
          IF notlast THEN
            BEGIN
              p:=p+z*h[ij+n];
              h[ij+n]:=h[ij+n]-p*r
            END;
          h[ij]:=h[ij]-p*q;
          h[ij-n]:=h[ij-n]-p
        END;{i}
      { accumulate transformations }
      FOR i:=low TO upp DO
        BEGIN
          ij:=i+n*k;
          p:=x*vecs[ij-n]+y*vecs[ij];
          IF notlast THEN
            BEGIN
              p:=p+z*vecs[ij+n]; vecs[ij+n]:=vecs[ij+n]-p*r
            END;
          vecs[ij]:=vecs[ij]-p*q;
          vecs[ij-n]:=vecs[ij-n]-p
        END;{i}

cont3:
    END; { k }
  GOTO nextit;
  { one root found }

onew:
  wr[en]:=x+t; h[en+n*Pred(en)]:=wr[en]; wi[en]:=0.0;
  cnt[en]:=its; en:=na; GOTO nextw;
  { two roots found }

twow:
  p:=(y-x)/2.0; q:=sqr(p)+w; z:=ABS(q); z:=Sqrt(z);
  x:=x+t; h[en+n*Pred(en)]:=x; h[na+n*Pred(na)]:=y+t;
  cnt[en]:=-its; cnt[na]:=its;
  IF q > 0.0 THEN
    BEGIN  { real pair }
      IF p < 0 THEN z:=p-z ELSE z:=p+z;
      wr[na]:=x+z; s:=x-w/z; wr[en]:=s; wi[na]:=0.0; wi[en]:=0.0;
      x:=h[en+n*Pred(na)]; r:=Sqrt(sqr(x)+sqr(z)); p:=x/r; q:=z/r;
      FOR j:=na TO n DO
        BEGIN
          ij:=n*Pred(j);
          z:=h[na+ij]; h[na+ij]:=q*z+p*h[en+ij];
          h[en+ij]:=q*h[en+ij]-p*z
        END; { row modification }
      FOR i:=1 TO en DO
        BEGIN
          ij:=i+n*Pred(na); z:=h[ij]; h[ij]:=q*z+p*h[i+n*Pred(en)];
          ij:=i+n*Pred(en); h[ij]:=q*h[ij]-p*z
        END; { column modification }
      FOR i:=low TO upp DO
        BEGIN
          ij:=i+n*Pred(na); z:=vecs[ij];
          vecs[ij]:=q*z+p*vecs[i+n*Pred(en)];
          ij:=i+n*Pred(en); vecs[ij]:=q*vecs[ij]-p*z
        END { accumulate }
    END  { pair of real roots }
             ELSE
    BEGIN  { complex pair }
      wr[na]:=x+p; wr[en]:=wr[na]; wi[na]:=z; wi[en]:=-z;
    END; { two roots found }
  en:=en-2; GOTO nextw;

  { all roots found, now backsubstitute }
fin:
IF NOT Schur THEN {Eigenvektoren vollstaendig}
  BEGIN
  norm:=0.0; k:=1;
  FOR i:=1 TO n DO
    BEGIN
      FOR j:=k TO n DO norm:=norm+ABS(h[i+n*Pred(j)]);
      k:=i
    END; { norm }
  { backsubstitution }
  FOR en:=n DOWNTO 1 DO
    BEGIN
      p:=wr[en]; q:=wi[en]; na:=Pred(en);
      IF q = 0.0 THEN
        BEGIN  { real vektor }
          m:=en; h[en+n*Pred(en)]:=1.0;
          FOR i:=na DOWNTO 1 DO
            BEGIN
              ij:=n*Pred(en);
              w:=h[i+n*Pred(i)]-p; r:=h[i+ij];
              FOR j:=m TO na DO r:=r+h[i+n*Pred(j)]*h[j+ij];
              IF wi[i] < 0.0 THEN
                BEGIN z:=w; s:=r END
                             ELSE
                BEGIN
                  m:=i;
                  IF wi[i] = 0.0 THEN
                    IF w <> 0.0 THEN h[i+ij]:=-r/w
                                ELSE h[i+ij]:=-r/(macheps*norm)
                                 ELSE
                    BEGIN  { Loesung des folgenden Gleichungssystems }
                    {        [ w  x ] [ h[i,en]   ]   [ -r ]      }
                    {  solve |      |*|           | = |    |      }
                    {        [ y  z ] [ h[i+1,en] ]   [ -s ]      }
                      ij:=i+n*Pred(i);
                      x:=h[ij+n];y:=h[Succ(ij)];
                      q:=Sqr(wr[i]-p)+Sqr(wi[i]);
                      ij:=i+n*Pred(en);
                      t:=(x*s-z*r)/q; h[ij]:=t;
                      IF ABS(x)>ABS(z)THEN h[Succ(ij)]:=(-r-w*t)/x
                                      ELSE h[Succ(ij)]:=(-s-y*t)/z
                    END {wi[i]>0.0}
                END{wi[i]>=0.0}
            END{i}
        END{real vector}
                ELSE
        IF q<0.0 THEN
          BEGIN {complex vector associated with lambda=p-i*q}
            m:=na;
            IF ABS(h[en+n*Pred(na)])>ABS(h[na+n*Pred(en)])
                   THEN
              BEGIN
                h[na+n*Pred(na)]:=-(h[en+n*Pred(en)]-p)/
                                             h[en+n*Pred(na)];
                h[na+n*Pred(en)]:=-q/h[en+n*Pred(na)]
              END
                   ELSE CDIV(-h[na+n*(en-1)],0.0,h[na+n*(na-1)]-p,q,
                         h[na+n*(na-1)],h[na+n*(en-1)]);
            h[en+n*Pred(na)]:=1.0; h[en+n*Pred(en)]:=0.0;
            FOR i:=Pred(na) DOWNTO 1 DO
              BEGIN {i}
                w:=h[i+n*Pred(i)]-p; ra:=h[i+n*Pred(en)]; sa:=0.0;
                FOR j:=m TO na DO
                  BEGIN
                    ij:=i+n*Pred(j);
                    ra:=ra+h[ij]*h[j+n*Pred(na)];
                    sa:=sa+h[ij]*h[j+n*Pred(en)]
                  END;
                IF wi[i] < 0.0 THEN
                  BEGIN z:=w; r:=ra; s:=sa END
                               ELSE
                  BEGIN
                    m:=i;
                    IF wi[i]=0.0 THEN CDIV(-ra,-sa,w,q,
                                          h[i+n*(na-1)],h[i+n*(en-1)])
                                 ELSE
                      BEGIN {Loesung der komplexen Gleichungen}
         { [w+q*i    x]  [ h[i,na]+h[i,en]*i        ]   [ -ra-sa*i  ]}
         { |          |  |                          | = |           |}
         { [  y  z+q*i]  [ h[i+1,na]+h[i;1,en] *i   ]   [  -r-s*i   ]}
                        ij:=i+n*Pred(i);
                        x:=h[ij+n]; y:=h[Succ(ij)];
                        vr:=Sqr(wr[i]-p)+Sqr(wi[i])-q*q;
                        vi:=(wr[i]-p)*2.0*q;
                        IF (vr=0.0) AND (vi=0.0) THEN
                          vr:=macheps*norm*(ABS(w)+ABS(q)+ABS(x)
                              +ABS(y)+ABS(z));
                        CDIV(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi,
                             h[i+n*(na-1)],h[i+n*(en-1)]);
                        IF ABS(x)>ABS(z)+ABS(q) THEN
                          BEGIN
                            ij:=i+n*Pred(na);
                            h[Succ(ij)]:=(-ra-w*h[ij]+
                                           q*h[i+n*Pred(en)])/x;
                            ij:=i+n*Pred(en);
                            h[Succ(ij)]:=(-sa-w*h[ij]-
                                           q*h[i+n*Pred(na)])/x
                          END
                                                ELSE
                          BEGIN
                            ij:=i+n*Pred(na);
                            CDIV(-r-y*h[ij],-s-y*h[i+n*(en-1)],z,q,
                                 h[ij+1],h[i+1+n*(en-1)])
                          END
                      END {wi[i]>0.0}
                  END {wi[i]>=0.0}
              END {i}
          END {complex vector}
    END {backsubstitution};
  {vectors of isolated roots}
  FOR i:=1 TO Pred(low) DO
    FOR j:=Succ(i) TO n DO vecs[i+n*Pred(j)]:=h[i+n*Pred(j)];
  FOR i:=Succ(upp) TO n DO
    FOR j:=Succ(i) TO n DO vecs[i+n*Pred(j)]:=h[i+n*Pred(j)];
  {multiply by transformation matrix to give vectors of original
   full matrix}
  FOR j:=n DOWNTO low DO
    BEGIN
      IF j<=upp THEN m:=j ELSE m:=upp;l:=Pred(j);
      IF wi[j] < 0.0 THEN
        BEGIN
          FOR i:=low to upp DO
            BEGIN
              y:=0.0;z:=y;
              FOR k:=low TO m DO
                BEGIN
                  ij:=i+n*Pred(k);y:=y+vecs[ij]*h[k+n*Pred(l)];
                  z:=z+vecs[ij]*h[k+n*Pred(j)]
                END;
              vecs[i+n*Pred(l)]:=y;vecs[i+n*Pred(j)]:=z
            END {i}
        END
                   ELSE
        IF wi[j] = 0.0 THEN
          FOR i:=low TO upp DO
            BEGIN
              z:=0.0;
              FOR k:=low TO m DO
                z:=z+vecs[i+n*Pred(k)]*h[k+n*Pred(j)];
              vecs[i+n*Pred(j)]:=z
            END {i}
    END {j};
  END {Eigenvektoren vollstaendig}
                  ELSE
  BEGIN {Schur-Vektoren}
    j:=1; { Nullsetzen unterhalb Schurform }
    REPEAT
      IF wi[j] = 0 THEN
        FOR i:=Succ(j) TO n DO h[i+n*Pred(j)]:=0.0
                   ELSE
        BEGIN
          FOR i:=j+2 TO n DO
             BEGIN h[i+n*Pred(j)]:=0.0; h[i+n*j]:=0.0 END;
          j:=Succ(j)
        END;
      j:=Succ(j);
    UNTIL j >= n
  END; {Schur-Vektoren}
fail:
END; { hqr2s }

{--------------------------------------------------------------------}


PROCEDURE EINORM (n: INTEGER; VAR Vec_in, Vec_out, wi_in);


VAR
  vecin:  RealMat Absolute Vec_in;
  vecout: RealMat Absolute Vec_out;
  wi:     RealVek Absolute wi_in;
  i,j,k,l,ij: INTEGER;
  betrag,maximum,maxre,maxim: Extended;

BEGIN
  j:=1; { j-ter Eigenwert, d.h. j-te Spalte der Eigenvektormatrix }
  REPEAT
    ij:=Pred(j)*n;
    IF wi[j] = 0.0 THEN { Normierung eines reellen Eigenvektors }
      BEGIN
        maximum:=0.0; k:=1;
        FOR i:=1 TO n DO { Suche des betragsgroeszten Elementes }
          BEGIN
            betrag:=ABS(vecin[i+ij]);
            IF betrag > maximum THEN
              BEGIN
                k:=i; { Merken der aktuellen Zeilen-Nr. }
                maximum:=betrag { neues Maximum }
              END
          END; { i }
        maximum:=vecin[k+ij];  { vorzeichenrichtiges Vektorelement }
                               { mit dem groeszten Betrag          }
        IF maximum <> 0.0 THEN { Division zur Normierung }
          FOR i:=1 TO n DO
            vecout[i+ij]:=vecin[i+ij]/maximum
                          ELSE { Nullvektor }
          FOR i:=1 TO n DO
            vecout[i+ij]:=0.0;
        Inc(j)  { Naechste Spalte der Eigenvektormatrix }
      END
                 ELSE { Normierung eines komplexen Eigenvektors }
      BEGIN
        maximum:=0.0; k:=1;
        FOR i:=1 TO n DO { Suche des beragsgroeszten Elementes }
          BEGIN
            betrag:=CABS(vecin[i+ij],vecin[i+ij+n]);
            IF betrag > maximum THEN
              BEGIN
                k:=i; { Merken der aktuellen Zeilen-Nr. }
                maximum:=betrag { neues Maximum }
              END
          END;  { i }
        maxre:=vecin[k+ij];    { vorzeichenrichtiger Realteil des }
                               { betragsgroeszten Vektorelementes }
        maxim:=vecin[k+ij+n];  { zugehoeriger Imaginaerteil       }
        IF maximum <> 0.0 THEN { komplexe Div. zur Normierung }
          FOR i:=1 TO n DO
            CDIV (vecin[i+ij],vecin[i+ij+n],maxre,maxim,vecout[i+ij],
                                                      vecout[i+ij+n])
                          ELSE { komplexer Nullvektor }
          FOR l:=j TO Succ(j) DO { Real- und Imaginaerteil in vecout }
            FOR i:=1 TO n DO   { alle n Zeilen des Vektors }
              vecout[i+Pred(l)*n]:=0.0;
        Inc(j,2)   { Nachste Spalte der Eigenvektormatrix }
      END;
  UNTIL j > n   { alle Spalten der Eigenvektormatrix fertig }
END; { EINORM }

{--------------------------------------------------------------------}


PROCEDURE EIGRF (N: INTEGER; VAR A_io,V_out,wr_out,wi_out);


VAR
   A: RealMat Absolute A_io;
   V: RealMat Absolute V_out;
   wr: RealVek Absolute wr_out;
   wi: RealVek Absolute wi_out;
   bal,cnt: ^RealVek;
   low,upp: INTEGER;
   vecsize: LongInt;

BEGIN
(*  if severe_error then Exit;  *)
  vecsize:=n*RealSize;
(*  if maxavail < 2*vecsize then
    begin
     GL_Fehler(9,severe_error); Exit
    end;*)
  getmem(bal,vecsize);
  getmem(cnt,vecsize);
  BALANCE(N,A,low,upp,bal^);
  ORTHES(N,low,upp,A,wr);
  ORTRANS(N,low,upp,A,V,wr);
  HQR2S(N,low,upp,A,V,wr,wi,cnt^,FALSE);
(*  if not severe_error then *)
    begin
      BALBAK(N,low,upp,N,bal^,V); EINORM(N,V,V,wi)
    end;
  freemem(cnt,vecsize);
  freemem(bal,vecsize)
END; (* EIGRF *)

{--------------------------------------------------------------------}

PROCEDURE LIAPS (n: INTEGER; VAR A_in,Q_io,P_out);


LABEL  kstep,one,two,cont1,cont2;
TYPE   Reroot = ARRAY[1..MaxVekSize] OF BOOLEAN;
VAR
   A:  RealMat Absolute A_in;
   Q:  RealMat Absolute Q_io;
   P:  RealMat Absolute P_out;
   AQ: ^RealMat;
   v1,v2: ^RealVek;
   wi: Reroot;
   i,j,k,l,f,g,ij: INTEGER;
   matsize,vecsize: LongInt;
   x,y,z,Akk: Extended;
   Korr,Bast: BOOLEAN;

PROCEDURE PSOLVE (VAR A,AQ: RealMat; VAR v2: RealVek;
                  VAR wi: Reroot; VAR n,k: INTEGER;
                                  VAR x,y,z:  Extended);

VAR  i,l,ij: INTEGER; sl,slp1,d,d1,d2,d3,d4: Extended;

BEGIN { PSOLVE }
  l:=Pred(k);
  REPEAT
    Inc(l); sl:=v2[l]; ij:=n*Pred(l);
    FOR i:=1 TO Pred(l) DO
      sl:=sl-(y*A[i+ij]+z*AQ[i+ij])*v2[i];
    IF wi[l] THEN
      v2[l]:=sl/(x+y*A[l+ij]+z*AQ[l+ij])
                      ELSE
      BEGIN
        slp1:=v2[Succ(l)];
        FOR i:=1 TO Pred(l) DO
          slp1:=slp1-(y*A[i+ij+n]+z*AQ[i+ij+n])*v2[i];
        d1:=x+y*A[l+ij]+z*AQ[l+ij];
        d2:=y*A[Succ(l)+ij]+z*AQ[Succ(l)+ij];
        d3:=y*A[l+ij+n]+z*AQ[l+ij+n];
        d4:=x+y*A[Succ(l)+ij+n]+z*AQ[Succ(l)+ij+n];
        d:=1.0/(d1*d4-d2*d3);
        v2[l]:=d*(d4*sl-d2*slp1); l:=Succ(l);
        v2[l]:=d*(d1*slp1-d3*sl)
      END; {if}
  UNTIL l >= n;
END; { PSOLVE }


BEGIN { LIAPS }
 (* if severe_error then Exit; *)
  matsize:=sqr(n)*RealSize;
  vecsize:=n*RealSize;
(*  if maxavail < (matsize+2*vecsize) then
    begin
      GL_Fehler(11,severe_error); Exit
    end; *)
  getmem(AQ,matsize);

  getmem(v1,vecsize);
  getmem(v2,vecsize);
  Korr:=true; i:=0; wi[n]:=TRUE;
  REPEAT
    Inc(i);
    IF ABS(A[Succ(i)+n*Pred(i)])=0 THEN wi[i]:=TRUE
                                   ELSE
      BEGIN wi[i]:=FALSE; Inc(i); wi[i]:=FALSE END;
  UNTIL i>=Pred(n);
  FOR i:=1 TO n DO
    BEGIN
      f:=Pred(i); IF f=0 THEN f:=1;
      FOR j:=f TO n DO
        BEGIN
          g:=Succ(j); IF g>n THEN g:=n; x:=0.0;
          FOR l:=f TO g DO x:=x+A[i+n*Pred(l)]*A[l+n*Pred(j)];
          AQ^[i+n*Pred(j)]:=x
        END { j }
    END; { i }
  k:=0; Bast:=TRUE;

kstep:
  Inc(k); IF k>n THEN GOTO cont1; Akk:=A[k+n*Pred(k)];
  IF wi[k] THEN GOTO one ELSE GOTO two;
one:
  FOR i:=1 TO Pred(k) DO v2^[i]:=P[k+n*Pred(i)];
  FOR i:=k TO n DO
    BEGIN
      x:=-Q[i+n*Pred(k)]; l:=Pred(k); ij:=n*Pred(k);
      FOR j:=1 TO l DO x:=x-A[j+ij]*P[i+n*Pred(j)];
      v2^[i]:=x
    END; { i }
  x:=Akk; y:=1.0; z:=0.0;

  PSOLVE(A,AQ^,v2^,wi,n,k,x,y,z);

  IF Bast THEN
    FOR i:=1 TO n DO P[i+n*Pred(k)]:=v2^[i]
          ELSE
    FOR i:=1 TO n DO
      BEGIN ij:=n*Pred(k); Q[i+ij]:=P[i+ij]; P[i+ij]:=v2^[i] END;
  GOTO kstep;
two:
  FOR i:=1 TO Pred(k) DO
    BEGIN
      x:=Q[i+n*Pred(k)]; l:=Pred(k);
      FOR j:=1 TO l DO x:=x+A[j+n*Pred(k)]*P[i+n*Pred(j)];
      v1^[i]:=x
    END; { i }
  z:=A[Succ(k)+n*Pred(k)];
  FOR i:=k TO n DO
    BEGIN
      x:=Q[i+n*Pred(k)]; y:=Q[i+n*k]; l:=Pred(k); ij:=n*Pred(k);
      FOR j:=1 TO l DO
        BEGIN
          x:=x+A[j+ij]*P[i+n*Pred(j)]; y:=y+A[j+ij+n]*P[i+n*Pred(j)]
        END; { j }
      v1^[i]:=x; v2^[i]:=z*y
    END; { i }
  FOR i:=1 TO Pred(k) DO v2^[i]:=P[k+n*Pred(i)];
  y:=A[Succ(k)+n*k];
  FOR i:=k TO Pred(n) DO
    BEGIN
      x:=v2^[i]; ij:=n*Pred(i);
      FOR j:=1 TO Succ(i) DO x:=x-A[j+ij]*v1^[j];
      v2^[i]:=x-y*v1^[i]
    END; { i }
  x:=v2^[n];
  FOR j:=1 TO n DO x:=x-A[j+n*Pred(n)]*v1^[j];
  v2^[n]:=x-y*v1^[n];
  x:=Akk*y-A[k+n*k]*A[Succ(k)+n*Pred(k)]; y:=Akk+y; z:=1.0;

  PSOLVE(A,AQ^,v2^,wi,n,k,x,y,z);

  y:=-1.0/A[Succ(k)+n*Pred(k)];
  FOR i:=1 TO Pred(k) DO v1^[i]:=P[Succ(k)+n*Pred(i)];
  FOR i:=k TO Pred(n) DO
    BEGIN
      x:=v1^[i]; ij:=n*Pred(i);
      FOR j:=1 TO Succ(i) DO x:=x+A[j+ij]*v2^[j];
      v1^[i]:=y*(x+Akk*v2^[i])
    END; { i }
  x:=v1^[n]; ij:=n*Pred(n);
  FOR j:=1 TO n DO x:=x+A[j+ij]*v2^[j];
  v1^[n]:=y*(x+Akk*v2^[n]);
  l:=Succ(k);
  IF Bast THEN
    FOR i:=1 TO n DO
      BEGIN P[i+n*Pred(k)]:=v2^[i]; P[i+n*Pred(l)]:=v1^[i] END
          ELSE
    FOR i:=1 TO n DO
      BEGIN
        ij:=i+n*Pred(k); Q[ij]:=P[ij]; P[ij]:=v2^[i];
        ij:=i+n*Pred(l); Q[ij]:=P[ij]; P[ij]:=v1^[i]
      END; {i,if}
  k:=l; GOTO kstep;

cont1:
  IF NOT Korr THEN GOTO cont2;
  IF NOT Bast THEN
    BEGIN
      FOR i:=1 TO n DO
        BEGIN
          ij:=i+n*Pred(i); P[ij]:=P[ij]+Q[ij];
          FOR j:=Succ(i) TO n DO
            BEGIN
              ij:=i+n*Pred(j);
              P[j+n*Pred(i)]:=P[ij]+Q[ij]; P[ij]:=P[ij]+Q[ij]
            END { j }
        END; { i }
      GOTO cont2;
    END; {if}
  FOR i:=1 TO n DO
    FOR j:=1 TO n DO
      BEGIN
        x:=0.0; l:=Succ(j); IF j=n THEN l:=n;
        FOR k:=1 TO l DO x:=x+P[i+n*Pred(k)]*A[k+n*Pred(j)];
        IF i=j THEN
          BEGIN ij:=i+n*Pred(i); Q[ij]:=Q[ij]+2.0*x END
               ELSE
          BEGIN
            ij:=i+n*Pred(j); Q[j+n*Pred(i)]:=Q[ij]+x; Q[ij]:=Q[ij]+x
          END
      END; {j,i}
  Bast:=FALSE; k:=0; GOTO kstep;
cont2:
  freemem(AQ,matsize);
  freemem(v1,vecsize);
  freemem(v2,vecsize);
END; { LIAPS }

{--------------------------------------------------------------------}

PROCEDURE UTCU (n: INTEGER; VAR U_in,C_io; rueck: BOOLEAN);

{   UTCU berechnet C*=U'CU mit U' als transponiertem U und C als
    symmetrischer Matrix. C ist zerlegbar in C=T+T' mit T als oberer
    Dreiecksmatrix. Damit laeszt sich C* aus C*=U'TU+(U'TU)' berech-
    nen, wobei wegen der Dreiecksgestalt von T gegenueber der direk-
    ten Ermittlung von C*=U'CU Multiplikationen eingespart werden.  }
VAR
  U: RealMat Absolute U_in;
  C: RealMat Absolute C_io;
  h: ^RealVek;
  i,j,k,ij,ji: INTEGER;
  s: Extended;
  vecsize: LongInt;

BEGIN (* UTCU *)
(*  if severe_error then Exit; *)
  vecsize:=n*RealSize;
(*  if maxavail < vecsize then
    begin
      GL_Fehler(10,severe_error); Exit
    end; *)
  getmem(h,vecsize);
  IF rueck THEN TRANSP(n,U);
  FOR i:=1 TO n DO   { wegen C=T+T' Diagonalelemente halbieren }
    BEGIN ij:=i+n*Pred(i); C[ij]:=0.5*C[ij] END;
  FOR i:=1 TO n DO  { Multiplikation T U }
    BEGIN
      FOR j:=1 TO n DO
        BEGIN
          s:=0.0;
          FOR k:=i TO n DO s:=s+C[i+n*Pred(k)]*U[k+n*Pred(j)];
          h^[j]:=s
        END;  { Zeile der Ergebnismatrix fertig }
      FOR j:=1 TO n DO C[i+n*Pred(j)]:=h^[j]
    END; { Matrixprodukt }
  FOR j:=1 TO n DO  { Multiplikation U'(TU) }
    BEGIN
      FOR i:=1 TO n DO
        BEGIN
          s:=0.0;
          FOR k:=1 TO n DO s:=s+U[k+n*Pred(i)]*C[k+n*Pred(j)];
          h^[i]:=s
        END;  { Spalte der Ergebnismatrix fertig }
      FOR i:=1 TO n DO C[i+n*Pred(j)]:=h^[i]
    END;  { Matrixprodukt }
  FOR i:=1 TO n DO  { C*=U'TU+(U'TU)' }
    FOR j:=i TO n DO
      BEGIN
        ij:=i+n*Pred(j); ji:=j+n*Pred(i);
        C[ij]:=C[ij]+C[ji]; C[ji]:=C[ij]
      END;
  IF rueck THEN TRANSP(n,U);
  freemem(h,vecsize)
END;  { UTCU }

(* ---------------------------------------------------------------- *)

PROCEDURE ROCOEX (n: INTEGER; VAR A_io);

{ Fuer die quadratische Matrix A wird eine Orthogonaltransformation
  derart vorgenommen, dasz die Zeilen und die Spalten jeweils in um-
  gekehrter Reihenfolge der Elemente in A stehen. Das entspricht der
  Vor- u. Nachmultiplikation von A mit einer Transformationsmatrix T,
  deren Nebendiagonalelemente 1 und alle uebrigen Elemente 0 sind.  }

VAR
  A: RealMat Absolute A_io;
  ii,jj: INTEGER;
  s: REAL;

BEGIN
  ii:=1; jj:=sqr(n);
  REPEAT
    s:=A[ii]; A[ii]:=A[jj]; A[jj]:=s; Inc(ii); Dec(jj)
  UNTIL ii >= jj
END; { ROCOEX }


(* ---------------------------------------------------------------- *)
(* 5.1.2. Procedure-Bibliothek der Basismodule ------------------- *)

PROCEDURE Ljapunow (N: INTEGER; VAR A_io,Q_io,P_out,U_io; LjpForm2:BOOLEAN;VAR skip_io: BOOLEAN);

LABEL fail;
(*skip_Io=TRUE: A hat bereits Schurform und Schurvektoren U liegen vor
  skip_iO=False:A ist instabil
    LjpForm2=True : AP + PA' = -Q
    LjpForm2=FALSE: A'P + PA = -Q *)
  Type RealVek = Array [1..20] of real;
VAR
  A: Mat Absolute A_io;
  Q: Mat Absolute Q_io;
  P: Mat Absolute P_out;
  U: Mat Absolute U_io;
Skip: BOOLEAN Absolute Skip_IO;
  wr,wi,cnt: ^RealVek;
  j,low,upp: INTEGER;
  vecsize: LongInt;

BEGIN  (* Ljapunow *)
 // if severe_error then Exit;
  vecsize:=n*RealSize;
  getmem(wr,vecsize); getmem(wi,vecsize);
  getmem(cnt,vecsize);
  IF NOT skip THEN  {zunaechst Schur-Trans.}
    BEGIN
      low:=1; upp:=N;
      ORTHES(N,low,upp,A,wr^);
      ORTRANS(N,low,upp,A,U,wr^);
      HQR2S(N,low,upp,A,U,wr^,wi^,cnt^,TRUE);
      Skip:=True;{Stabilitaetstest}

      if Not Skip then Goto fail;
    //  if severe_error then Goto fail;
    END; { if }
  IF LjpForm2 THEN { Aschur:= T A'schur T }
  BEGIN TRANSP(N,A); ROCOEX(N,A) END;
  UTCU(N,U,Q,FALSE); IF LjpForm2 THEN ROCOEX(N,Q);
  LIAPS(N,A,Q,P); IF LjpForm2 THEN ROCOEX(N,P);
  UTCU(N,U,P,TRUE);
fail:
  freemem(cnt,vecsize);
  freemem(wi,vecsize);
  freemem(wr,vecsize)
END;  (* Ljapunow *)


PROCEDURE Lia_Probe(n: INTEGER; VAR A_in,Q_io,P_in; LjpForm2: BOOLEAN);

VAR
     A: RealMat Absolute A_in;
     Q: RealMat Absolute Q_io;
     P: RealMat Absolute P_in;
     i,j,k: INTEGER; s: Extended;
BEGIN
  IF LjpForm2 THEN TRANSP(n,A);
  FOR i:=1 TO n DO
    FOR j:=1 TO n DO
      BEGIN  { Liapunow-Gleichung mit eingesetzter Loesung }
        s:=Q[i+n*Pred(j)];
        FOR k:=1 TO n DO
          s:=s+A[k+n*Pred(i)]*P[k+n*Pred(j)]
              +P[i+n*Pred(k)]*A[k+n*Pred(j)];
        Q[i+n*Pred(j)]:=s
      END;
  IF LjpForm2 THEN TRANSP(n,A)
END; { Lia_Probe }

{* ---------------------------------------------------------------- *}

Procedure Eigenwerte (N:Integer;Var Ag_io;Var Ew_io:Vek;Var TM_io; Var Stabil_io:Boolean);
  Type RealVek = Array [1..20] of real;
Var
  Ag: RealMat Absolute Ag_io;
  TM: RealMat Absolute TM_io;
  EW: Vek     Absolute EW_io;
  Stabil: Boolean   Absolute Stabil_io;

  wr,wi,cnt: ^RealVek;
  j,low,upp: INTEGER;
  vecsize: LongInt;
  ETm,ED : Real;
Begin
With ARBMAT^ do begin
  vecsize:=N*RealSize;
  getmem(wr,vecsize);
  getmem(wi,vecsize);
  getmem(cnt,vecsize);
  low:=1; upp:=N;
  ORTHES(N,low,upp,Ag,wr^);

  ORTRANS(N,low,upp,Ag,TM,wr^);
  HQR2S(N,low,upp,Ag,TM,wr^,wi^,cnt^,TRUE);

  Stabil:=True;{Stabilitaetstest}
 (* ETm:=0;ED:=0;
  For j:=low To upp Do begin
      ETm:=ETm+1/Abs(wr^[j]);ED:=ED+1/Sqrt(SQR(wr^[j])+SQR(wi^[j]));
      EW[j+1]:=wr^[j];EW[N+2+j]:=wi^[j];
      If wr^[j]>=0 Then Stabil:=False;
  End;
      EW[1]:=ETm;EW[N+2]:=ETm/ED; *)
  freemem(cnt,vecsize);
  freemem(wi,vecsize);
  freemem(wr,vecsize);
End;
End;



End.
