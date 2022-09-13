unit uDaten;

{$mode ObjFPC}{$H+}

interface

 uses
  Classes, SysUtils;

CONST
MaxMatSize = 2500;
MaxVekSize = 50;
RealSize = 8;
VekDim=16;
MatDim =1000;
DUH =49;   {U,H Dimension 7x7 muss Hesse-Dimension einschliessen!}
EPS = 1E-6; {Abruchschranke in VLAK,REGULA }

TYPE
RealMat  = ARRAY[1..MaxMatSize] OF REAL;
RealVek  = ARRAY[1..MaxVekSize] OF REAL;
IntVek   = ARRAY[1..MaxVekSize] of Integer;
Vek  = ARRAY[1..VekDim] OF REAL;
Mat  = ARRAY[1..MatDim] OF REAL;
MUH  = ARRAY[1..DUH] OF REAL;

HeapOutMat = RECORD
   Ag       : Mat;         {Systemmatrix des geschlossenen Kreises}
   PX,PY    : Mat;         {Loesung der Ljapunow-Gleichungen}
   QX,QY    : Mat;         {Parametermatrizen dito}

 piT,giT,gi,pi: Mat;       {BFGS, DFP}
END;

HeapArbMat = RECORD
ZM,TM,SU  : MAT; {Hilfsmatrizen}
END;


VAR

{Zustandsraum-Regler-Matrizen}
          Ar: Mat; {Nr,Nr Systemmatrix}
          Br: Mat; {Nr,Ms Eingangsmatrix}
          Cr: Mat; {Ps,Nr Ausgangsmatrix}
          Dr: Mat; {Ms,Ps stat.Rückführmatrix}
          Rr: Mat; {Nr,Nr dyn.Kovarianzmatrix}
          Qr: Mat; {Nr,Nr dyn.Zustandswichtungsmatrix}
{einheitl.Optimierungsmatrizen im Berechnungsalgorithmus}
           A: Mat; { N, N Systemmatrix}
           B: Mat; { N, M Eingangsmatrix}
           C: Mat; { P, N Ausgangsmatrix}
           D: Mat; { P, M Durchgangsmatrix fuer Matlab}
           K: Mat; { M, P Reglermatrix => Nk=m*p beachte Hesse-Dimension!}
           Q: Mat; { N, N Zustangswichtungsmatrix}
{Hilfsmatrizen}
          Qx: Mat; { N, N Zustangswichtungsmatrix}
          Qz: Mat; { N, N Zustangswichtungsmatrix}
          Px: Mat; { N, N Ljapunov loesung}
          Py: Mat; { N, N Ljapunov loesung}
          Pz: Mat; { N, N Ljapunov loesung}
          Ax: Mat; { N, N Ljapunov loesung}
          Az: Mat; { N, N Ljapunov loesung}
          Ag: Mat; { N, N Geschlossener Regelkreis}
          TM: Mat; { N, N Eigenwerte}
          HM: Mat; { N, N Eigenwerte}
          HZ: Mat; { N, N Eigenwerte}
          Mt: Mat; { N, N Transponierte Matrix}
           R: Mat; { M, M Stellenergiewichtungsmatrix}
           V: Mat; { N, N Kovarianzmatrix x0}
           W: Mat; { P, P Kovarianzmatrix y0}
           U: Mat; { N, N Hesse-Transform }
          SU: Mat; { N, N Schur U-Matrix}
          He: Mat; { N, N Hesse-Transform }
           H: MUH; {Nk,Nk Hesse-Matrix Achtung DIM-Explosionsgefahr!}
     //      U: MUH; {Nk,Nk Gütefuntions-Wichtungsmatrix}


 Ky_LOG     : Mat; {freie Reglerelemente = 0 ; fixierten <> 0}

     //  InpMat :^HeapInpMat;
       OutMat :^HeapOutMat;
       ArbMat :^HEAPArbMat;

{Definition von Vektorgroessen Dimension:=Anzahl der freien Variablen}

     Ki_Succ: Mat;
         K_y: Mat;
      Gi_Vek: Mat;
     Gi_Pred: Mat;
      hi_Vek: Mat; {Nur in Gradient verwaltet}
          Ew: Vek;
          Pm: Vek;
          Wr: Vek;
          Wi: Vek;
{Definition von Skalaren-Variablen}

N,M,P      : Integer;     {Zustandsraum-Dimensionen}
Ns,Ms,Ps   : Integer;
Nk,Nr      : Integer;     {Anzahl der freien Parameter des Reglers}
VecSize    : LongInt;
LjpIter    : Integer;
LakIter    : Integer;
LjpForm2,skip   : Boolean;
Stabil     : Boolean;
MaInvSel   : Integer=2;
QDif       : Real;
low,upp    : Integer;

implementation

//{$I VLAKBIB.PAS}

end.
