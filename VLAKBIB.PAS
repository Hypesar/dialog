(* 5.1.2. Procedure-Bibliothek der Basismodule ------------------- *)

Procedure InHyp (QM : Mat;ZP,SQ :Integer;
             Var ZM : Mat;ZN,SM :Integer;
                          Z0,S0 :Integer);
Var k,l,a0,b0 : Integer;
Begin
{$R+}
  a0:=Z0+ZN*(S0-1)-1;
  For l:=1 to SQ do Begin
      b0:=ZP*(l-1);
      For k:=1 to ZP do ZM[a0+k]:=QM[b0+k];
      a0:=a0+ZN;
  End;
{$R-}
End;

Procedure AusHyp (QM : Mat;ZP,SQ :Integer;
              Var ZM : Mat;ZN,SM :Integer;
                           Z0,S0 :Integer);
Var k,l,a0,b0 : Integer;
Begin
{$R+}
  a0:=Z0+ZP*(S0-1)-1;
  For l:=1 to SQ do Begin
      b0:=ZN*(l-1);
      For k:=1 to ZN do ZM[b0+k]:=QM[a0+k];
      a0:=a0+ZP;
  End;
{$R-}
End;

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


Procedure MASUB (VAR A_i,B_i,C_o:MAT;N,M:Integer);
var j:Integer;
Begin
  For j:=1 to N*M do C_O[j]:=A_I[j]-B_I[j];
End;

Procedure MASUH (VAR A_i:MUH;Var B_i:Mat;Var C_o:MUH;N,M:Integer);
var j:Integer;
Begin
  For j:=1 to N*M do C_O[j]:=A_I[j]-B_I[j];
End;

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
       For i:=1 TO N*N Do QX[i]:=Q[i]+QX[i];

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
      For i:=1 TO sqr(n) Do QY[i]:=V[i]+QY[i];
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

Procedure A_BKC (N,M,P: Integer;VAR A,B,K,C,Ag: MAT);
Var  i: Integer;
Begin
  MMAMA(B,K,Ag,n,m,p);
  MMAMA(Ag,C,Ag,n,p,n);
  For i:=1 TO sqr(n) Do Ag[i]:=A[i]-Ag[i];
{  Maus2d (Ag,N,N,'Ag');}
End; { A_BKC }

Function Funktional : Extended;
Begin
With ArbMat^,OutMat^ do begin
  Vek_zu_Mat (K_y,K,Ky_Log,NK,M,P);
  QY_Aus (V,B,K,W,TM,N,M,P);
  MMAMA (PX,TM,HM,N,N,N);
{ k'Uk Modifizierung }

  TM:=K_y;
  MTRSP (TM,Nk,1);
  MMAMA (TM,U,ZM,1,Nk,Nk);
  MMAMA (ZM,K_y,TM,1,Nk,1);

  Funktional:=Spur(Hm,N)+TM[1];
END;
END;

Procedure Phi_x ;
VAR  J: Integer;
     Fx,Fy,Fu :Extended;
    SU,AG,QX,QY,Px : Mat;
   FORM2,Skip     : BOOLEAN;
BEGIN
  With InpMAT^,ArbMat^ do begin
    Vek_zu_Mat (K_y,K,Ky_Log,NK,M,P);
    A_BKC (N,M,P,A,B,K,C,Ag);
    QX:=Q;
    Skip:=False;{Kaltberechnung von P}
    Form2:=False;{A`P+PA=-Q}
    Ljapunow(N,Ag,QX,PX,SU,Form2,skip);
    MMAMA (PX,V,HM,N,N,N);
    Fx:=0;for j:=1 to N*N do fx:=fx+hm[j];
(*
    if skip Then Begin
      QY_Aus (V,B,K,W,QY,N,M,P);
      Skip:=True;{Kaltberechnung von P}
      Form2:=True;{AP+P'A=-Q}
      Ljapunow(N,Ag,QY,PX,SU,Form2,skip);
    End;
    if Skip Then Begin gotoxy(50,2);Write ('Matrix Ag :    Stabil');end
            Else Begin gotoxy(50,2);Write ('Matrix Ag :  Instabil');end;
*)
     gotoxy(50,5);write ('Phi x     : ',Fx:12);

  end;
end;

Procedure Stab_Set;

Begin
  If SW_Stab=0 Then With ArbMat^,OutMat^,InpMat^ do begin
  (* A-Stabilisierung setzen *)
    SW_Stab:=1;
    Einheits_Matrix (Tm,N);
    MAFAK (SW_Stab,Tm,N,N);
    MASUB (A,Tm,A,N,N);
    write (chr(7));
  End;
End;

(*--------------- Problemspezifische Proceduren -----------------------*)

(* 5.1.3. *)

PROCEDURE LJAPXY (Var PX,PY : MAT;VAR A_STABIL : BOOLEAN);
VAR  J: Integer;
   FORM2,Skip     : BOOLEAN;
BEGIN
  With InpMAT^,ArbMat^,OutMat^ do begin
    QX_Aus (Q,C,K,R,QX,N,M,P);
    A_BKC (N,M,P,A,B,K,C,Ag);
    Skip:=False;{Kaltberechnung von P}
    Form2:=False;{A`P+PA=-Q}
    Ljapunow(N,Ag,QX,PX,SU,Form2,skip);
     A_Stabil:=Skip;
(*
    if Skip Then Begin gotoxy(50,2);Write ('Matrix Ag :    Stabil');end
            Else Begin gotoxy(50,2);Write ('Matrix Ag :  Instabil');end;
    gotoxy(50,5);write ('LJAP-Iter : ',LjpIter:8);
*)
    If Skip Then Begin
      QY_Aus (V,B,K,W,QY,N,M,P);
      Skip:=True;{Kaltberechnung von P}
      Form2:=True;{AP+P'A=-Q}
      Ljapunow(N,Ag,QY,PY,SU,Form2,skip);
    End;
     LjpIter:=LjpIter+1;
  end;
end;

(* 5.1.6. *)

Procedure GRADIENT (Alf_i:Extended;Var Ki_i,Ki_Succ,Gi_o:Mat;
                               Var Skalar_io:Extended);
  var gi,hi         : Extended;
      J             : Integer;
Procedure GRAD_KUHN (Var Ki_i,Gi_o:Mat);
{Problemspezifischer Gradient}
Var G,Z1 : Mat;
Begin
With INPMAT^,OUTMAT^,ARBMAT^ do begin
  Vek_zu_Mat (Ki_i,K,Ky_Log,NK,M,P);
{ -2B`XYC`}
  TM:=B;
  MTRSP (TM,N,M);
  MMAMA (TM,PX,HM,M,N,N);
  MMAMA (HM,PY,TM,M,N,N);
  HM:=C;
  MTRSP (HM,P,N);
  MMAMA (TM,HM,G,M,N,P);
  MAFAK (-2,G,M,P);
{ 2B`XBKW}
  TM:=B;
  MTRSP (TM,N,M);
  MMAMA (TM,PX,HM,M,N,N);
  MMAMA (HM,B,TM,M,N,M);
  MMAMA (TM,K,HM,M,M,P);
  MMAMA (HM,W,TM,M,P,P);
  MAFAK (-2,TM,M,P);
  MASUB (G,TM,G,M,P);
{ 2RKCYC`}
  MMAMA (R,K,HM,M,M,P);
  MMAMA (HM,C,TM,M,P,N);
  MMAMA (TM,PY,HM,M,N,N);
  TM:=C;
  MTRSP (TM,P,N);
  MMAMA (HM,TM,ZM,M,N,P);
  MAFAK (-2,ZM,M,P);
  MASUB (G,ZM,G,M,P);

  Mat_zu_Vek (G,Gi_o,Ky_Log,NK,M,P);
{ 2Uk (Modifizierung)}

  MMAMA (U,Ki_i,ZM,Nk,Nk,1);
  MAFAK (-2,ZM,Nk,1);
  MASUB (Gi_o,ZM,Gi_o,Nk,1);

End; {with}
End; {KUHN}

Begin
With ArbMat^,InpMat^ do begin
  If Alf_i=0 Then begin
     Grad_KUHN (Ki_i,Gi_o);
     MMAMA (H,Gi_o,hi_Vek,Nk,Nk,1);
     Ki_Succ:=Ki_i;
    end else begin
     ZM:=hi_Vek;
     MAFAK (Alf_i,ZM,Nk,1);
     MASUB (Ki_i,ZM,Ki_Succ,Nk,1);
     Grad_KUHN (Ki_Succ,Gi_o);
  End;
  Skalar_io:=0;
  For j:=1 to Nk do Skalar_io:=Skalar_io-Gi_o[j]*hi_Vek[j];
END;
End;

(* 5.1.6. *)

Procedure BFGS_Formel (Var H : MUH;Var gi_i,pi_i:Mat;Nk: Integer);
var piTgi,giTHigi : Extended;
    x,y           : Extended;
    j             : Integer;
    piT,giT       : Mat;
{Hi+1=Hi-(Higipi'+pigi'Hi)/(pi'gi)+pigi'[1+gi'Higi/(pi'gi)]/(pi'gi)}
Begin
  With ArbMat^ do begin
    piT:=pi_i;giT:=gi_i;
    MTRSP (piT,Nk,1);MTRSP (giT,Nk,1);
{ pi'gi }
    MMAMA (piT,gi_i,ZM,1,Nk,1);
    piTgi:=ZM[1];
    y:=1/piTgi;
{ gi'Higi }
    MMAMA (giT,H,ZM,1,NK,NK);
    MMAMA (ZM,gi_i,TM,1,Nk,1);
    giTHigi:=TM[1];
    x:=giTHigi/piTgi;
{ Hgipi' * Y }
    MMAMA (H,gi_i,TM,Nk,Nk,1);
    MMAMA (TM,piT,ZM,Nk,1,Nk);
    MAFAK (y,ZM,Nk,Nk);
    MASUH (H,ZM,H,Nk,Nk);
{ pigi'Hi * Y }
    MMAMA (pi_i,giT,TM,Nk,1,Nk);
    MMAMA (TM,H,ZM,Nk,Nk,Nk);
    MAFAK (y,ZM,Nk,Nk);
    MASUH (H,ZM,H,Nk,Nk);
{ pipi' *  -Y }
    MMAMA (pi_i,piT,ZM,Nk,1,Nk);
    MAFAK (-y,ZM,Nk,Nk);
    MASUH (H,ZM,H,Nk,Nk);
    MAFAK (x,ZM,Nk,Nk);
    MASUH (H,ZM,H,Nk,Nk);
    LakIter:=lakIter+1;
  END;
End;

(* 5.1.6. *)

Procedure DFP_Formel (Var H:MUH;Var gi_i,pi_i:Mat;Nk: Integer);
var piTgi,giTHigi : Extended;
    x,y           : Extended;
    j             : Integer;
    piT,giT       : Mat;
{ Hi+1 := Hi - Higigi'Hi/gi'Higi + pipi'/pi'gi }
Begin
  With ArbMat^ do begin
    piT:=pi_i;giT:=gi_i;
    MTRSP (piT,Nk,1);MTRSP (giT,Nk,1);
{ gi'Higi }
    MMAMA (giT,H,TM,1,Nk,Nk);
    MMAMA (TM,gi_i,ZM,1,Nk,1);
    giTHigi:=ZM[1];
{ Higigi'Hi }
    MMAMA (H,gi_i,ZM,Nk,Nk,1);
    MMAMA (ZM,TM,HM,Nk,1,NK);
    MAFAK (1/giTHigi,HM,Nk,Nk);
    MASUH (H,HM,H,Nk,Nk);
{ pi'gi }
    MMAMA (piT,gi_i,ZM,1,Nk,1);
    piTgi:=ZM[1];
{ pipi' }
    MMAMA (pi_i,piT,ZM,Nk,1,NK);
    MAFAK (-1/piTgi,ZM,Nk,Nk);
    MASUH (H,ZM,H,Nk,Nk);
    LakIter:=lakIter+1;
  END;
End;

(* 5.1.7. *)

Procedure KMIN  (Var Alfa:Extended;Var Ki_i,Ki_Succ,Gi_Vek:Mat;
                        VAR Grad:Extended);

Function Rk_instabil (Var Alfa,Grad:Extended):Boolean;
Var
           SU,Ag : Mat;
  FORM2,Skip     : BOOLEAN;
Begin
  With InpMAT^,OutMat^,ArbMat^ do begin
{ K:=Ki_Succ(Alfa) }
  Gradient (Alfa,K_y,Ki_Succ,Gi_Vek,Grad);
  VEK_zu_MAT (Ki_Succ,K,Ky_log,Nk,M,P);
  LJAPXY (PX,PY,A_Stabil);
  Rk_instabil:= not A_Stabil;
  Gradient (Alfa,K_y,Ki_Succ,Gi_Vek,Grad);
  End;{With}
End;

Procedure Fkt_Grad (Var Alfa,Fkt,Grad:Extended);
Var fai:Extended;
Form2,Skip:Boolean;
Begin
  With ArbMat^,OutMat^,InpMat^ do begin
{ K:=Ki_Succ(Alfa) }
    Gradient (Alfa,K_y,Ki_Succ,Gi_Vek,Grad);
    VEK_zu_MAT (Ki_Succ,K,Ky_log,Nk,M,P);
    LJAPXY (PX,PY,A_Stabil);
    Gradient (Alfa,K_y,Ki_Succ,Gi_Vek,Grad);
    QY_Aus (V,B,K,W,TM,N,M,P);
    MMAMA (PX,TM,HM,N,N,N);
{ k'Uk Modifizierung }
    TM:=Ki_Succ;
    MTRSP (TM,Nk,1);
    MMAMA (TM,U,ZM,1,Nk,Nk);
    MMAMA (ZM,Ki_Succ,TM,1,Nk,1);
    fai:=Spur(Hm,N)+TM[1];
    Fkt:=fai;
  END;
END;

Function Restart (Var Gi_Vek:Mat;Var H : MUH):Extended;
Begin
  With ArbMat^,OutMat^,InpMat^ do begin
    TM:=Gi_Vek;
    MAFAK (-1,TM,Nk,Nk);
    MMAMA (H,TM,HM,Nk,Nk,1);
    TM:=Gi_Vek;
    MTRSP (TM,Nk,1);
    MMAMA (TM,HM,ZM,1,Nk,1);
    Restart:=ZM[1];
  End;
END;

Procedure Stab_Kor;
  Var SW_END : Extended;
Begin
  With ArbMat^,OutMat^,InpMat^ do begin
  If SW_Stab<>0 Then Begin
(* A putzen *)
    Einheits_Matrix (Tm,N);
    MAFAK (-SW_Stab,Tm,N,N);
    MASUB (A,Tm,A,N,N);
    If SW_Stab<1E-3 Then Sw_Stab:=0 Else Begin
      OpenWindow(22,10,60,15);
      HeadLine(12,'Stabilisierung');
          str(SW_Stab:12,Alfs0);
          Kette36:=Concat('SW_STAB :',alfs0);
          printat (2,1,2,Kette36);
      Einheits_Matrix (Tm,N);
      SW_Stab:=SW_Stab/10;
      MAFAK (+SW_Stab,Tm,N,N);
      MASUB (A,Tm,A,N,N);
      VEK_zu_MAT        (K_y,K,Ky_log,Nk,M,P);
      A_BKC (N,M,P,A,B,K,C,Ag);
      Eigenwerte (N,Ag,Ew,A_Stabil);
      SW_END:=SW_STAB*10;
      If Not A_Stabil Then
        While NOT A_Stabil do Begin
          Einheits_Matrix (Tm,N);
          MAFAK (-SW_Stab,Tm,N,N);
          MASUB (A,Tm,A,N,N);
          SW_Stab:=(SW_END+SW_Stab)/2;
          Einheits_Matrix (Tm,N);
          MAFAK (+SW_Stab,Tm,N,N);
          MASUB (A,Tm,A,N,N);
          A_BKC (N,M,P,A,B,K,C,Ag);
          Eigenwerte (N,Ag,Ew,A_Stabil);
          str(SW_Stab:12,Alfs0);
          Kette36:=Concat('SW_STAB :',alfs0);
          printat (2,3,2,Kette36);
        End Else Begin
{         Maus2D (Ew,N+1,2,'Ew (A-BKC)');}
{
         While A_Stabil and (SW_Stab>1e-3) do Begin
          Einheits_Matrix (Tm,N);
          MAFAK (-SW_Stab,Tm,N,N);
          MASUB (A,Tm,A,N,N);
          Einheits_Matrix (Tm,N);
          SW_Stab:=SW_Stab/2;
          MAFAK (+SW_Stab,Tm,N,N);
          MASUB (A,Tm,A,N,N);
          A_BKC (N,M,P,A,B,K,C,Ag);
          Eigenwerte (N,Ag,Ew,A_Stabil);
          str(SW_Stab:12,Alfs0);
          Kette36:=Concat('SW_STAB :',alfs0);
          printat (2,3,2,Kette36);
         End;
         If NOT A_Stabil Then Begin
          Einheits_Matrix (Tm,N);
          MAFAK (-SW_Stab,Tm,N,N);
          MASUB (A,Tm,A,N,N);
          SW_Stab:=SW_Stab*2;
          Einheits_Matrix (Tm,N);
          MAFAK (+SW_Stab,Tm,N,N);
          MASUB (A,Tm,A,N,N);
          A_BKC (N,M,P,A,B,K,C,Ag);
          Eigenwerte (N,Ag,Ew,A_Stabil);
         End;
 }
       End;

      (* In jedem Stabilisierungsschritt Restart *)
      Einheits_Matrix(H,Nk);
      Alfa:=0;
      Gradient (Alfa,K_y,Ki_Succ,Gi_Vek,Grad);

      CloseWindow;
    End;
 End;
 End;
End;

Label
      sch100,
      sch200,
      sch300,
      sch400;

var ai,ai0,ai1,ai2,ai3,aih,hi1,hi2,gi1,gi2,gih,fi1,fi2,fih:Extended;
    j : Integer;
    alfs0,alfs1,alfs2,alfs3: String[12];
    Kette34 : String[22];
    Kette36 : String[34];
Begin

  With ArbMat^,OutMat^,InpMat^ do begin
(* DFP-Quasi-Newton-Verfahren nach Kuhn Bild 4.1 *)
    Alfa:=0;
    Gradient (Alfa,K_y,Ki_Succ,Gi_Vek,Grad);

    aih:=1;
    Fi_Old:=0;
    OpenWindow(2,4,39,23);
    HeadLine(10,'Optimierungsalgorithmus');
    FootLine(10,'ESC = Unterbrechung'); 
   REPEAT
    SCINIT;
    Gi_Pred:=Gi_Vek;
    Dfr:=Abs((Fi_New-Fi_Old)/Fi_New);
    Fi_Old:=Fi_New;

(* Beginn der Eindimensionalen Minimumsuche *)
(*
    If (Dfr<1E-5) and (Auswahl2=1) Then begin
         Alfa:=1;
         Gradient (Alfa,K_y,Ki_Succ,Gi_Vek,Grad);
         VEK_zu_MAT (Ki_Succ,K,Ky_log,Nk,M,P);
         A_BKC(N,M,P,A,B,K,C,Ag);
         Eigenwerte (N,Ag,Ew,A_stabil);
         if A_stabil then goto sch400;
    end;
*)
    Alfa:=aih;
    printat(2,3,2,'SCHRITT 1 :           Gradient    ');

    If Alfa<=0 Then Alfa:=0.5;
    ai0:=2*Alfa;ai1:=0;ai3:=0;
Sch100:
    if RK_instabil (ai0,gi2) then begin
      ai3:=ai0;
      str(ai3:12,Alfs3);str(gi2:12,alfs0);
      Kette36:=Concat('Alpha 3 :',alfs3,' ',alfs0);
      printat (2,6,2,Kette36);	
      
      if ai1=0.1*ai0 then goto Sch200 else begin
        ai0:=0.1*ai0;
        goto Sch100;
      end;
     end else if gi2>0 then begin
               ai2:=ai0;
               str(ai2:12,Alfs2);str(gi2:12,alfs0);
               Kette36:=Concat('Alpha 2 :',alfs2,' ',alfs0);
               printat (2,5,2,Kette36);
               goto Sch300;{*}
               end else begin
                ai1:=ai0;
                str(ai1:12,Alfs1);str(gi2:12,alfs0);
                Kette36:=Concat('Alpha 1 :',alfs1,' ',alfs0);
               printat (2,4,2,Kette36);*)
                gi1:=gi2;
                if ai3=10*ai0 then goto Sch200 else begin
                  ai0:=10*ai0;
                  goto Sch100;
                end;
              end;

Sch200:
    printat(2,8,2,'SCHRITT 2 :           Gradient    ');
    ai2:=Sqrt(ai1*ai3);
    if RK_instabil (ai2,gi2) then begin
      str(ai2:12,Alfs2);str(gi2:12,alfs0);
      Kette36:=Concat('Alpha 2 :',alfs2,' ',alfs0);
      printat (2,10,2,Kette36);
      ai3:=ai2;
      str(ai3:12,Alfs3);str(gi2:12,alfs0);
      Kette36:=Concat('Alpha 3 :',alfs3,' ',alfs0);
      printat (2,11,2,Kette36);
      goto Sch200;
     end else if gi2>0 then Begin
                str(ai2:12,Alfs2);str(gi2:12,alfs0);
                Kette36:=Concat('Alpha 2 :',alfs2,' ',alfs0);
                printat (2,10,2,Kette36);
                goto Sch300;
               End else begin
                ai1:=ai2;
                gi1:=gi2;
                str(ai1:12,Alfs1);str(gi2:12,alfs0);
                Kette36:=Concat('Alpha 1 :',alfs1,' ',alfs0);
                printat (2,9,2,Kette36);
                str(ai2:12,Alfs2);str(gi2:12,alfs0);
                Kette36:=Concat('Alpha 2 :',alfs2,' ',alfs0);
                printat (2,10,2,Kette36);
                goto Sch200;
              end;
sch300:
      FKT_GRAD (ai1,fi1,gi1);
      FKT_GRAD (ai2,fi2,gi2);
      printat(2,14,2,'SCHRITT 3 : Hermite Interpolation ');
      printat(2,15,2,'          Index 1       Index 2   ');
      REPEAT
      ai0:=aih;
    (* Hermite Interpolation *)
       hi1:=3*(fi1-fi2)/(ai2-ai1)+gi1+gi2;
       hi2:=(hi1*hi1-gi1*gi2);
       If hi2>0 Then Begin
         hi2:=SQRT(hi2);
         aih:=ai2-(ai2-ai1)*(hi1+gi2-hi2)/(2*hi1+gi1+gi2);
         FKT_GRAD (aih,fih,gih);
         If gih>0 THEN Begin
             ai2:=aih;gi2:=gih;fi2:=fih;
          end ELSE Begin
             ai1:=aih;gi1:=gih;fi1:=fih;
         End;
        End Else Begin
        aih:=(ai1+ai2)/2;ai0:=aih;
       End;
       str(ai1:12,Alfs1);str(ai2:12,Alfs2);
       Kette36:=Concat('Alpha  :',alfs1,'  ',alfs2);
       printat (2,16,2,Kette36);
       str(gi1:12,Alfs1);str(gi2:12,Alfs2);
       Kette36:=Concat('Gradient',alfs1,'  ',alfs2);
       printat (2,17,2,Kette36);
       str(fi1:12,Alfs1);str(fi2:12,Alfs2);
       Kette36:=Concat('Funktion',alfs1,'  ',alfs2);
       printat (2,18,2,Kette36); 
     Until (ABS((aih-ai0)/(aih+ai0))<0.05) or (Abs(ai1-ai2)<Arfeps);{End Hermite}

   (* Quasi -Newton *)
sch400:
       Alfa:=aih;
       Gradient (Alfa,K_y,Ki_Succ,Gi_Vek,Grad);

       DKyn:=0;Dkyz:=0;
       For j:=1 to Nk do Begin
         Dkyn:=Dkyn+SQR(K_y[j]-Ki_Succ[j]);
         Dkyz:=Dkyz+SQR(Ki_Succ[j]);
       End;
       Dky:=Sqrt(Dkyn/Dkyz);
       For j:=1 to Nk do Begin
         gi[j]:=Gi_Vek[j]-Gi_Pred[j];{gi:=gi2-gi1}
         pi[j]:=Ki_Succ[j]-K_y[j];   {pi:=ki2-ki1}
       End;
       K_y:=Ki_Succ;
       DFP_Formel  (H,gi,pi,Nk);
       If Restart (Gi_Vek,H)>0 Then begin
         Einheits_Matrix(H,Nk);
         Alfa:=0;
         Gradient (Alfa,K_y,Ki_Succ,Gi_Vek,Grad);
       End;
       Stab_KOR;
(*       ChangeWindow (2); Entwurfsinformationen *)
       SCINIT;
       SWPOLE:=1;

       Alfa:=0;
       FKT_GRAD (Alfa,fi_New,Grad);
       SWPOLE:=0;
       Gi_Pred:=Gi_Vek;
       str(fih:12,Alfs0);str(grad:12,Alfs1);
       Kette36:=Concat('Fkt:',alfs0,' Grad:',Alfs1);
       printat (2,3,2,Kette36);
       str(LakIter:5,Alfs0);str(LjpIter:6,Alfs1);
       Kette36:=Concat('BFGS-Iter:',alfs0,'  ','Ljap.-Ber.:',Alfs1);
       printat (2,4,2,Kette36);

(*       ChangeWindow (3);  Optimierungsinformationen *)
       Read;
       if Keypressed Then ESC_CODE:=ORD(READKEY);
   Until (Abs(Grad)<Arfeps) or (ESC_CODE=esc_s) ;
(*   Until (Abs(Grad)<Arfeps) *)
   CloseWindow; 
  End;{With}
END;{KMIN}

(* 5.1.8. *)

Procedure Init_Mat;

Begin
 With ARBMAT^ do begin
  ArfEps:=Eps;
  Einselemente_Matrix (Ky_Log,N*N); {alle Kij freie Parameter}
  K_y:=K;Nr:=0;
  Nk:=M*P;
  Einheits_Matrix(Qs,Ns);
  Einheits_Matrix(Rs,Ms);
  Null_Matrix(U,Nk*Nk);
  Einheits_Matrix(H,Nk);
  Einheits_Matrix(V,Ns);
  Null_Matrix(W,Ps*Ps);
  Q:=Qs;R:=Rs;

 End;
End;

Procedure Start;
Begin
  With InpMAT^,OutMat^,ArbMat^ do begin
  EINLSD (3,3,Ns,Ms,Ps,As,Bs,Cs,Dr,' Zustandsraumsystem ',
                  '---> Matrix As','---> Matrix Bs',
                  '---> Matrix Cs','---> Matrix Ds');
     Mbez:='K';MATIN(4,13,Ms,Ps,Dr,'Startmatrix K ');
  Ag:=As;
  A:=As;B:=Bs;C:=Cs;K:=Dr;
  N:=Ns;M:=Ms;P:=Ps;Nr:=0;
  Nk:=M*P;
  SW_Stab:=0;
  INIT_MAT; {Wichtungsmatrizen Q,R,V,W,U = I ;Kij->K_y; Nk=M*P }
  END;
  Auswahl1:=2;
END;

(* 5.1.9. *)

Procedure Statisch(Var Auswahl:Integer);
CONST
  M2: ARRAY[1..5] Of SelektString =
     (' Neustart ZR-Mat.:A,B,C,D&K einl. ',
      ' Reglerstartmatrix K einlesen     ',
      ' Integrale Streckenweiterung      ',
      ' Stabilisierung A:=A-I     mit K=0',
      ' Riccati-Regler-Entwurf C=I setzen');

Var  E:Mat;I:Byte;
begin
  With InpMAT^,OutMat^,ArbMat^ do begin
  Auswahl:=1;
  OpenWindow(20,12,57,21);
  Headline(8,' M2-Startoperationen ');
  Footline(16,' ESC => M1->2: ');
  Selekt(2,2,M2,5,Auswahl);
  REPEAT
    Case Auswahl of
    1: Begin
         Start;
         Auswahl2:=2;
         Auswahl:=esc_s
       End;
    2: if Nr=0 Then Begin
         Mbez:='K';MATIN(4,13,Ms,Ps,Dr,'--->K (Dr)');
         K:=Dr;
         K_y:=K;
         Auswahl:=esc_s
       End Else begin
         write (chr(7));
         Auswahl2:=1;
         Auswahl:=esc_s;
       End;
    3: If (Nr=0) and (N<>0) Then Begin
         N:=Ns+Ps;P:=Ps+Ps;M:=Ms;
         Null_Matrix  (A,N*N);
         InHyp        (As,Ns,Ns,A,N,N,1,1);
         InHyp        (Cs,Ps,Ns,A,N,N,Ns+1,1);
         Null_Matrix  (B,N*M);
         InHyp        (Bs,Ns,Ms,B,N,M,1,1);
         Null_Matrix  (C,P*N);
         InHyp        (Cs,Ps,Ns,C,P,N,1,1);
         InHyp        (Cs,Ps,Ns,C,P,N,Ps+1,1);
         Einheits_Matrix (E,Ps);
         InHyp        (E,Ps,Ps,C,P,N,Ps+1,Ns+1);
         Null_Matrix  (K,M*P);
{        InHyp        (Dr,Ms,Ps,K,M,P,1,1);}
(*Streckenmatrizen vom ZR-ArbMat=>ZRs-Mat zurÅck bringen *)
         As:=A;Bs:=B;Cs:=C;Dr:=K;
         Ns:=N;Ps:=P;Ms:=M;
         Init_Mat; {Wichtungsmat. etc. mit neuen Dimensionen }
         Stab_Set;
         Auswahl:=esc_s
       End;
    4: Begin
         Stab_Set;
         Auswahl:=esc_s
       End;
    5: Begin
(*Streckenmatrizen vom ZR-ArbMat=>ZRs-Mat zurÅck bringen *)
         N:=Ns;M:=Ms;P:=Ns;Ps:=Ns;
         Einheits_Matrix (C,N);
         Null_Matrix  (K,M*N);
         As:=A;Bs:=B;Cs:=C;Dr:=K;
         Init_Mat; {Wichtungsmat. etc. mit neuen Dimensionen }
         Auswahl:=esc_s
       End;
    6: ArfEps:=RREAD(20,18,'1E-15   ','1E-1    ',1E-15,1E-1);
    End;
  UNTIL Auswahl=esc_s;
{
  A_Stabil:=True;
  For i:=1 to M*P do If K[i]>0 Then A_Stabil:=False;
 If A_Stabil Then Begin
  Stab_Set;

  Einheits_Matrix (Tm,N);
  MASUB (A,Tm,A,N,N);
  A_BKC (N,M,P,A,B,K,C,Ag);
  Eigenwerte (N,Ag,Ew,A_Stabil);
  If A_Stabil Then Begin
    Einheits_Matrix (Tm,N);
    MAFAK (-1,Tm,N,N);
    MASUB (A,Tm,A,N,N);
    SW_Stab:=0;
    for i:=1 to N do If EW[i+1]>=SW_Stab Then SW_Stab:=EW[i+1];
    SW_Stab:=SW_Stab+2;
  (* A-Stabilisierung setzen *)
    Einheits_Matrix (Tm,N);
    MAFAK (SW_Stab,Tm,N,N);
    MASUB (A,Tm,A,N,N);
    write (chr(7));
  End;

 End;
}
 CloseWindow;
  End;{With}
END;

(* 5.1.10. *)

Procedure Abbruch_Epsilon;
CONST
  M4: ARRAY[1..1] Of SelektString =
     (' Abr.Eps 1E-6    ');
Var i,Auswahl : Integer;
BEGIN
With ARBMAT^ do begin
Auswahl:=1;
OpenWindow(20,12,70,15);
Headline(10,' M4-Abbruchepsilon ');
Footline(16,' ESC => M3->1: ');
REPEAT
Selekt(2,2,M4,1,Auswahl);
Case Auswahl of
  1 : ArfEps:=RREAD(20,18,'1E-20   ','1E-1    ',1E-20,1E-1);
END;
UNTIL AUSWAHL=esc_s;
CloseWindow;
End;
End;

(* 5.1.11. *)

Procedure Diagonal (Var Auswahl:Integer);
CONST
  M5: ARRAY[1..7] Of SelektString =
     (' Qs- Zustandswichtung     ',
      ' Rs- Stellenergiewichtung ',
      ' V - E{x`x} Anfangszustand',
      ' Ws- E{y`y} Me·stîrung    ',
      ' U - Konvergenz Wichtung  ',
      ' Qr- Zustandwichtung Regl.',
      ' Rr- Stellenergie Regler  ');
Var SYM:Boolean;
BEGIN
With ARBMAT^ do begin
Auswahl:=1;
OpenWindow(20,12,70,24);
Headline(8,' M5-Diagonalelemente ');
Footline(16,' ESC => M3->2: ');
SYM:=False;
REPEAT
Selekt(2,2,M5,7,Auswahl);
Case Auswahl of
1 : Begin
      Mbez:='Qs';MATIN(4,13,Ns,Ns,Qs,'--->Qs-Diag ');
      Q:=Qs;
      If Nr>0 Then begin
        Einheits_matrix (Q,N);
        InHyp (Qs,Ns,Ns,Q,N,N,1,1);
        InHyp (Qr,Nr,Nr,Q,N,N,Ns+1,Ns+1);
      End;
{      Maus2d (Q,N,N,'Q');}
    End;
2 : Begin
      Mbez:='Rs';MATIN(4,13,Ms,Ms,Rs,'--->Rs-Diag ');
      R:=Rs;
      If Nr>0 Then begin
        Einheits_Matrix (R,M);
        InHyp (Rs,Ms,Ms,R,M,M,1,1);
        InHyp (Rr,Nr,Nr,R,M,M,Ms+1,Ms+1);
      End;
{      Maus2d (R,M,M,'R');}
      End;

3 : Begin
      Mbez:='V';MATIN(4,13,N,N,V,'--->V-Diag ');
    End;
4 : Begin
      Einheits_Matrix(Ws,Ps);
      Mbez:='Ws';MATIN(4,13,Ps,Ps,Ws,'--->Ws-Diag ');
      W:=Ws;
      If Nr>0 Then begin
        Null_Matrix (W,P*P);
        InHyp (Ws,Ps,Ps,W,P,P,1,1);
      End;
{      Maus2d (W,P,P,'W');  }
    End;
5 : Begin
      Einheits_Matrix(U,Nk);
      Mbez:='U';MATIN(4,13,Nk,Nk,U,'--->U-Diag ');
{      Maus2d (U,Nk,Nk,'Ur');}
    End;
6 : If Nr>0 Then Begin
      Mbez:='Qr';MATIN(4,13,Nr,Nr,Qr,'--->Qr-Diag ');
      InHyp (Qr,Nr,Nr,Q,N,N,Ns+1,Ns+1);
{      Maus2d (Q,N,N,'Q');   }
    End;
7 : If Nr>0 Then Begin
      Mbez:='Rr';MATIN(4,13,Nr,Nr,Rr,'--->Rr-Diag ');
      InHyp (Rr,Nr,Nr,R,M,M,Ms+1,Ms+1);
{      Maus2d (R,M,M,'R');}
    End;
END;
UNTIL (AUSWAHL=esc_s) or NOT SYM;
CloseWindow;
End;
End;

(* 5.1.12. *)

Procedure Ky_Fixierung ;
Var i:Integer;
    E:Mat;
BEGIN
With ARBMAT^ do begin
OpenWindow(3,3,33,20);
Headline(3,' M6-Reglerelemente-Fixierung ');
Footline(10,' ESC => M3->3: ');
   Mbez:='Ky_Log';MATIN(4,13,Ms,Ps,E,'--->Ky_Log ');
if Nk<>0 Then begin
   InHyp (E,Ms,Ps,Ky_log,M,P,1,1);
End Else Ky_Log:=E;
MAT_zu_Vek (K,K_y,Ky_Log,Nk,M,P);
CloseWindow;
End;
End;

(* 5.1.13 *)

Procedure Vorgaben (Var Auswahl:Integer);
CONST
  M3: ARRAY[1..3] Of SelektString =
     (' M4-Abbruchepsilon                 ',
      ' M5-Diagonalelemente               ',
      ' M6-Reglerelementefixierung nur Dr');

BEGIN
With ARBMAT^ do begin
OpenWindow(20,12,58,19);
Headline(10,' M3-Entwurfsvariationen ');
Footline(16,' ESC => M1->3: ');
REPEAT
  Selekt(2,2,M3,3,Auswahl);
  Case Auswahl of
  1 : Abbruch_Epsilon;
  2 : Diagonal  (Auswahl33);
  3 : Ky_Fixierung ;
  End;
UNTIL auswahl=esc_s;
CloseWindow;
End;
END;

(* 5.1.14. *)

Procedure Variation;
BEGIN
  With InpMAT^,OutMat^,ArbMat^ do begin
  Auswahl3:=2;
  VORGABEN (Auswahl3);
  Vek_zu_Mat (K_y,K,Ky_Log,NK,M,P);
  LjpIter:=1;
  SCINIT;
  SWPOLE:=1;
  LJAPXY (PX,PY,A_Stabil);
  SWPOLE:=0;
  fi_New:=Funktional;
  fi_Min:=fi_New;
  Einheits_Matrix (H,Nk);{Fuer Gi_Skalar}
  Alf_Nul:=0;
  GRADIENT (Alf_Nul,K_y,Ki_Succ,Gi_Vek,Gi_Skalar);{OptimalKontrolle}
  str(fi_New:12,Alfs0);
  str(fi_New:12,Alfs0);str(Gi_Skalar:12,Alfs1);
  Kette36:=Concat('Fkt:',alfs0,' Grad:',Alfs1);
  printat (2,3,2,Kette36);
  str(LakIter:5,Alfs0);str(LjpIter:6,Alfs1);
  Kette36:=Concat('BFGS-Iter:',alfs0,'  ','Ljap.-Ber.:',Alfs1);
  printat (2,4,2,Kette36);
  END;
  Auswahl1:=3;
END;

(* 5.1.15. *)

Procedure Optimal;
var Alfa : Extended;
Begin
  With InpMAT^,OutMat^,ArbMat^ do begin

   Vek_zu_Mat (K_y,K,Ky_Log,NK,M,P);
   LJAPXY (PX,PY,A_Stabil);
   Einheits_Matrix (H,Nk);
   LakIter:=0;
   LjpIter:=0;
   Alf_Nul:=0;
   Alfa:=1;
   GRADIENT (Alf_Nul,K_y,Ki_Succ,Gi_Vek,Gi_Skalar);
   fi_New:=Funktional;
   QRU_Fak:=exp(Int(ln(fi_New)/ln(10))*ln(10));
{   MAFAK   (1/QRU_Fak,Q,N,N);
   MAFAK   (1/QRU_Fak,R,M,M);}
   Gi_Skalar:=1;
   Dky:=1;

   ESC_CODE:=0;
   KMIN (Alfa,K_y,Ki_Succ,Gi_Vek,Gi_Skalar);
   If Nr>0 Then Begin
      AusHyp (K,M,P,Dr,Ms,Ps,1,1);
     End Else Dr:=K;
   If Keypressed Then ESC_CODE:=ORD(READKEY);
   If ESC_Code=esc_s Then Auswahl1:=2 Else Auswahl1:=4;
  END;
{  MAFAK   (QRU_Fak,Q,N,N);
  MAFAK   (QRU_Fak,R,M,M);}
  Write (chr(7));
END;


(* 5.1.16 *)

Procedure Protokoll;
Type S20=String[20];
Procedure Transfer (Var A,B,C,D,K:Mat;N,M,P:Integer;DaNam:S20);

Var FVar : Text;
    Wort: ARRAY [1..9] of Char;
    SW  : String[8];
    I   : Integer;
  Stabil: Boolean;
Procedure D2M (Var M:Mat;Z,S:Integer;MaNam:Char);
Var    I,J,l :Integer;
       Amin,Amax,MFak,Fak:Real;
       SFak : String[6];
       Mi   :Mat;
Begin
    Mi:=M;
    Amin:=Mi[1];Amax:=Mi[1];
    For l:=2 to Z*S do Begin
       If Abs(Mi[l])<Abs(Amin) Then Amin:=Mi[l];
       If Abs(Mi[l])>Abs(Amax) Then Amax:=Mi[l];
    End;
    Mfak:=1;Fak:=1;
    If (Amin<>0)and(Amax<>0) Then Begin
     For l:=1 to 3 do Begin
      If (Abs(Amin)<Abs(1/Fak)) and (Abs(Amax)<Abs(1/Fak*1E3)) Then MFak:=Fak;
      Fak:=Fak*10;
     End;
    End;
    For l:=1 to 9 do Wort[l]:=' ';
    Wort[1]:=MaNam;
    Wort[2]:='=';
    Wort[3]:='[';
    If Mfak<>1 Then Begin
       For l:=1 to Z*S do Mi[l]:=Mi[l]*MFak;
       Str(1/MFak:5:3,SFak);
       SW:=ConCat (MaNam,'=',SFak);
       For l:=1 to 7 do Wort[l]:=Sw[l];
       Wort[8]:='*';
       Wort[9]:='[';
     End;
    Writeln(FVar,Wort);
    For i:=1 to z do Begin
       For j:=1 to s do Begin
         Wort[9]:=' ';
         Str(Mi[i+Z*(j-1)]:8:3,SW);
         For l:=1 to 8 do Wort[l]:=Sw[l];
         If j=s Then Begin
            If i=z Then Wort[9]:=']' Else Wort[9]:=';';
            Writeln(Fvar,Wort);
         End Else Write(Fvar,Wort);
       End;
     End;
     Writeln(Fvar);
End;
Begin
 With Arbmat^,OutMat^ do Begin
  Assign (FVar,DANAM);
  Rewrite(FVar);
  D2M (A,N,N,'A');
  D2M (B,N,M,'B');
  D2M (C,P,N,'C');
{  D2M (D,P,M,'D');}
  D2M (K,M,P,'K');
(* Diagonalelemente *)
  For  i:=1 to N do Tm[i]:=Q[i+N*(i-1)];
  D2M (Tm,1,N,'Q');
  For  i:=1 to M do Tm[i]:=R[i+M*(i-1)];
  D2M (Tm,1,M,'R');
  For  i:=1 to N do Tm[i]:=V[i+N*(i-1)];
  D2M (Tm,1,N,'V');
  if Nr>0 Then Begin
    Stabil:=False;
    For  i:=1 to Nk do Begin
       Tm[i]:=U[i+Nk*(i-1)];
       if Tm[i]<>0 Then Stabil:=True;
    End;
    if Stabil Then D2M (Tm,1,Nk,'U');
  End;
  A_BKC (N,M,P,A,B,K,C,Ag);
  Eigenwerte (N,Ag,Ew,Stabil);
  For i:=1 to 2*(N+1) do Tm[i]:=Ew[i];
  D2M (Tm,N+1,2,'E');
  Close(FVar);
 End;{With}
End;

CONST
  M7: ARRAY[1..5] Of SelektString =
     (' Ergebnistransfer =>KUARF.M f. MATLAB  ',
      ' Reglermatrix K [Dr Cr;-Br -Ar] Ky_log ',
      ' Ew (A-BKC)     Kriterien   Re=Tm Im=D ',
      ' Ew (Ar-BrDrCr) Kriterien   Re=Tm Im=D ',
      ' Arbeitsmatrizen A,B,C,K,R,Q,V,W und U ');
Var Auswahl7: Integer;
    Stabil  : Boolean;
Begin
With ARBMAT^ do begin
  Auswahl7:=1;
  OpenWindow(2,12,45,20);
  Headline(10,' M7-Protokollierung ');
  Footline(16,' ESC => M1->5: ');
  REPEAT
    Selekt(2,2,M7,5,Auswahl7);
    if Nk<>0 Then With InpMAT^,OutMat^,ArbMat^ do begin

      Case Auswahl7 of
      1 :  Begin
             Null_Matrix (D,P*M);
             Transfer (A,B,C,D,K,N,M,P,'\matlab\my\kuarf.m');
           End;
      2:   Begin
            Maus2D (K,M,P,'K');
            Maus2D (Ky_log,M,P,'Ky_log');

           End;
      3 :  Begin
           A_BKC (N,M,P,A,B,K,C,Ag);
           Eigenwerte (N,Ag,Ew,Stabil);
           Maus2D (Ew,N+1,2,'Ew (A-BKC)');
           End;
      4 :  If Nr<>0 Then Begin
           AusHyp (K,M,P,Dr,Ms,Ps,1,1);
           AusHyp (K,M,P,Cr,Ms,Nr,1,Ps+1);
           AusHyp (K,M,P,Br,Nr,Ps,Ms+1,1);
           AusHyp (K,M,P,Ar,Nr,Nr,Ms+1,Ps+1);
           MAFAK  (-1,Ar,Nr,Nr);
           MAFAK  (-1,Br,Nr,Ps);
           A_BKC  (Nr,Ps,Ms,Ar,Br,Dr,Cr,Ag);
           Eigenwerte (Nr,Ag,Ew,Stabil);
           Maus2D (Ew,Nr+1,2,'EW (Ar-BrDrCr)');
           Maus2D (Ar,Nr,Nr,'Ar');
           Maus2D (Br,Nr,Ps,'Br');
           Maus2D (Cr,Ms,Nr,'Cr');
           Maus2D (Dr,Ms,Ps,'Dr');
           End;
      5:   Begin
           Maus2D (A,N,N,'A');
           Maus2D (B,N,M,'B');
           Maus2D (C,P,N,'C');
           Maus2D (K,M,P,'K');
           Maus2d (Q,N,N,'Q');
           Maus2d (R,M,M,'R');
           Maus2d (V,N,N,'V');
           Maus2d (W,P,P,'W');
           Maus2d (U,Nk,Nk,'U');

           End;
      End;{case}
    End Else  Headline(7,' Noch Keine Matrizen da!');
  UNTIL auswahl7=esc_s;
  CloseWindow;
End;
END;

(* 5.1.17. *)

Procedure Dynamisch;
Var Arm,Arh,Crm,Crh,Nqi,l,i,j : Integer;
    E,Kfs : Mat;
Procedure I_Regler;
  var i : Integer;
  SWEW0 : Boolean;
Begin
 With InpMAT^,OutMat^,ArbMat^ do Begin
  If Pm[1]=1 Then Begin
    (* Asymmetrische Matrix Ar *)
     Arfeps:=1E-4;
     AusHyp  (Ky_Log,M,P,Kfs,Nr,1,Ms+1,P);
     For i:=1 to Nr do If PM[5+i]=0 Then Begin
        Kfs[i]:=0;
        Stab_Set;
     End;
     MAT_zu_Vek        (K,K_y,Ky_Log,Nk,M,P);
     InHyp   (Kfs,Nr,1,Ky_log,M,P,Ms+1,P);
    End Else Begin
    (* Symmetrische Matrix Ar *)
    If Nr=Ms Then Begin
      SWEW0:=True;
      For i:=1 to Nr do if PM[5+i]<>0 Then SWEW0:=False;
      If SWEW0 Then Begin
        Arfeps:=1E-6;
        Einheits_Matrix (E,Nr);         {Cr=i}
        InHyp           (E,Nr,Nr,K,M,P,1,Ps+1);
        Null_Matrix     (E,Nr*Nr);      {Ar=0}
        InHyp           (E,Nr,Nr,K,M,P,Ms+1,Ps+1);
        InHyp           (E,Nr,Nr,Ky_log,M,P,Ms+1,Ps+1);
        InHyp           (E,Nr,Nr,Ky_log,M,P,1,Ps+1);
        MAT_zu_Vek      (K,K_y,Ky_Log,Nk,M,P);
        Stab_Set;
      End;
    End;
  End;
 End;
End;

Begin
  With InpMAT^,OutMat^,ArbMat^ do begin
   OpenWindow(3,3,44,20);
   Headline(3,'M8-Reglerstrukturparameter');
   Footline(10,'   ESC=Ende ');
   printat(2,2,2,'Dynamische ARF : K=[Dr Cr;-Br -Ar]    ');
   printat(2,3,2,'symmetrisches-Ap[q ni;Ew]`asymmetrisch');
   printat(2,4,2,' ⁄ x 0 ø  x-freies Element    ⁄ 0 x ø ');
   printat(2,5,2,' ¿ 0 x Ÿ  Strukturparameter   ¿ 1 x Ÿ ');
   printat(2,6,2,'[2 1 1 0 0]   1. Spalte   [1 2 0 0 0] ');
   printat(2,7,2,'[-1 -1 0 0 0] 2. Spalte Startpole (=) ');
   printat(2,8,2,'Anzahl Elemente Nk=q*(m+p)+m*p <max 35');
   printat(2,9,2,'! q Û m bei sym. q-Ord. m-Inp. p-Outp.');
   Mbez:='PM ';MATIN(4,13,5,2,Pm,'Ap[q ni;Ew]`');
   CloseWindow;
   If Nr<>0 Then Begin
   (* Dr-Fixierung retten *)
      AusHyp            (Ky_log,M,P,Kfs,Ms,Ps,1,1);
      AusHyp            (K,M,P,Tm,Ms,Ps,1,1);
{      Maus2D (Kfs,ms,ps,'Kfs');}
      MAT_zu_Vek        (Tm,K_y,Ky_Log,Nk,M,P);
    End Else Begin
      Kfs:=Ky_log;
   End;

   Arfeps:=1E-3;
 if Round(PM[1]) in [1..Ms] Then Begin
   (*K,Ar,Cr,Ky_Log - Hyper-Matrizen kano.Form nach Kuhn aufbauen *)
   A:=As;B:=Bs;C:=Cs;K:=Dr;
   N:=Ns;M:=Ms;P:=Ps;{Stat. originaler Indizes Zustand der Matrizen}
   Nr:=0;
   For l:=1 to Round(PM[1]) do Nr:=Nr+Round(PM[l+1]);{Reglerordnung}
   Crh:=(M+Nr)*P;Arh:=Crh+M;{Adressen fÅr K,Ky_log}
   Crm:=0;Arm:=0;
   Nk:=Nr*(Ms+Ps)+Nk;   {Ms*Ps Anzahl der freien Kij}
   Null_Matrix (E,(Ns+Nr)*(Ns+Nr));
   A:=E;B:=E;C:=E;K:=E;
   Ar:=E;Br:=E;Cr:=E;
   Einselemente_Matrix  (Ky_log,(M+Nr)*(P+Nr));
   InHyp (Kfs,Ms,Ps,Ky_log,M+Nr,P+Nr,1,1);
   InHyp (E,Nr,Nr,Ky_log,M+Nr,P+Nr,Ms+1,Ps+1);
   Nqi:=0;
   For l:=1 To Round(PM[1]) do Begin {q-Anzahl Ari}
      (* Anzahl der Spalten von Aqi,Cqi *)
      For i:=1 To Round(PM[l+1]) do Begin
         (*Cr,Ky_log-Zeilen*)
         For j:=1 To M do Begin
            If l=j Then Begin
              Ky_log[Crh+j]:=0;
              If i=Round(PM[l+1]) Then Cr[Crm+j]:=1;
            End;
         End;
         Crm:=Crm+M;
         Crh:=Crh+M+Nr;
         (*Ar,Ky_log-Zeilen*)
         For j:=1 To Round(PM[l+1]) do Begin
            If j=i+1 Then Ar[Arm+j]:=1;
            (*letzte Spalte *)
            If i=Round(PM[l+1]) Then Begin
              Ky_log[Arh+Nqi+j]:=1;
                  Ar[Arm+Nqi+j]:=PM[5+Nqi+j];{Eigenwerte}
            End;
         End;
         Arm:=Arm+Nr;
         Arh:=Arh+M+Nr+Nqi;
      End;
      Nqi:=Nqi+Round(PM[l+1]);
    End;{l}
(* Indizes Transformation zu Hypermatrizen A,B,C,K *)
    N:=Ns+Nr;M:=Ms+Nr;P:=Ps+Nr;
(*Hyper-Matrizeninitialisierung*)
    Einheits_Matrix (E,Nr);
(* der Strecke *)
    InHyp           (As,Ns,Ns,A,N,N,1,1);
    InHyp           (Bs,Ns,Ms,B,N,M,1,1);
    InHyp           (E,Nr,Nr,B,N,N,Ns+1,Ms+1);
    InHyp           (Cs,Ps,Ns,C,P,N,1,1);
    InHyp           (E,Nr,Nr,C,P,N,Ps+1,Ns+1);
(* des Reglers *)
    InHyp           (Dr,Ms,Ps,K,M,P,1,1);
    MAFAK           (-1,Ar,Nr,Nr);
    InHyp           (Ar,Nr,Nr,K,M,P,Ms+1,Ps+1);
    InHyp           (Cr,Ms,Nr,K,M,P,1,Ps+1);
    I_Regler;
{   Maus2D (Ky_log,M,P,'Ky_log');}
(* der Wichtungsmatrizen *)
    Einheits_Matrix (Qr,Nr);Rr:=Qr;
    Null_Matrix     (R,M*M);
    InHyp           (Rs,Ms,Ms,R,M,M,1,1);
    InHyp           (Rr,Nr,Nr,R,M,M,Ms+1,Ms+1);
    Null_Matrix     (Q,N*N);
    InHyp           (Qs,Ns,Ns,Q,N,N,1,1);
    InHyp           (Qr,Nr,Nr,Q,N,N,Ns+1,Ns+1);
    Null_Matrix     (W,P*P);
    Einheits_Matrix (H,Nk);
    Null_Matrix     (U,Nk*Nk);
    Einheits_matrix (V,N);
    Auswahl1:=2;
   End Else Begin
     Write (chr(7));
     Auswahl1:=5;
  End;   {If}
  MAT_zu_Vek        (K,K_y,Ky_Log,Nk,M,P);

  LjpIter:=1;
  SCINIT;
  SWPOLE:=1;
  VEK_zu_MAT        (K_y,K,Ky_log,Nk,M,P);
  LJAPXY            (PX,PY,A_Stabil);
  SWPOLE:=0;
  Einheits_Matrix   (H,Nk);
  Alf_Nul:=0;
  GRADIENT          (Alf_Nul,K_y,Ki_Succ,Gi_Vek,Gi_Skalar);
  fi_New:=Funktional;
  fi_Min:=fi_New;
  str(fi_New:12,Alfs0);
  str(fi_New:12,Alfs0);str(Gi_Skalar:12,Alfs1);
  Kette36:=Concat('Fkt:',alfs0,' Grad:',Alfs1);
  printat (2,3,2,Kette36);
  str(LakIter:5,Alfs0);str(LjpIter:6,Alfs1);
  Kette36:=Concat('BFGS-Iter:',alfs0,'  ','Ljap.-Ber.:',Alfs1);
  printat (2,4,2,Kette36);
END;
END;
