unit urunmode;

{$mode ObjFPC}{$H+}

interface

uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, Menus, ExtCtrls,
  StdCtrls, Grids, uDaten, uLak, uTools, uViewInput;

type

  { TLjapunov }

  TLjapunov = class(TForm)
    btnCloseView: TButton;
    btnViewE: TRadioButton;
    btnQdif: TEdit;
    MainMenu1: TMainMenu;
    btnGrubLjaponv: TRadioGroup;
    BtnLjapunov1: TRadioButton;
    btnLjapunov2: TRadioButton;
    btnViewA: TRadioButton;
    btnViewAt: TRadioButton;
    btnViewP: TRadioButton;
    btnViewQ: TRadioButton;
    btnViewU: TRadioButton;
    btnViewSU: TRadioButton;
    btnViewEv: TRadioButton;
    btnJjapunovS2: TRadioButton;
    btnJjapunovS1: TRadioButton;
    btnMinv: TRadioButton;
    btnMaInv: TRadioButton;
    btnTmInv: TRadioButton;
    rbQrTrans: TRadioButton;
    rbEwHes: TRadioButton;
    rbQrthes: TRadioButton;
    rbBalance: TRadioButton;
    RadioGroup1: TRadioGroup;
    RadioGroup2: TRadioGroup;
    ViewGroup: TRadioGroup;
    StrGridLjapunov: TStringGrid;
    procedure btnCloseViewClick(Sender: TObject);
    procedure btnGrubLjaponvClick(Sender: TObject);
    procedure BtnLjapunov1Change(Sender: TObject);
    procedure btnLjapunov2Change(Sender: TObject);
    procedure btnMinvChange(Sender: TObject);
    procedure btnTmInvChange(Sender: TObject);
    procedure btnViewAChange(Sender: TObject);
    procedure btnViewAtChange(Sender: TObject);
    procedure btnViewEChange(Sender: TObject);
    procedure btnViewEvChange(Sender: TObject);
    procedure btnViewPChange(Sender: TObject);
    procedure btnViewQChange(Sender: TObject);
    procedure btnViewSUChange(Sender: TObject);
    procedure btnViewUChange(Sender: TObject);
    procedure btnQdifChange(Sender: TObject);
    procedure FormActivate(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure Label1Click(Sender: TObject);
    procedure btnJjapunovS2Change(Sender: TObject);
    procedure btnJjapunovS1Change(Sender: TObject);
    procedure btnMaInvChange(Sender: TObject);
    procedure rbBalanceChange(Sender: TObject);
    procedure rbEwHesChange(Sender: TObject);
    procedure rbQrthesChange(Sender: TObject);
    procedure rbQrTransChange(Sender: TObject);
  private

  public

  end;

var
  Ljapunov: TLjapunov;

implementation

{$R *.lfm}

{ TLjapunov }

procedure TLjapunov.FormActivate(Sender: TObject);
begin

end;

procedure TLjapunov.btnCloseViewClick(Sender: TObject);
begin
  Ljapunov.close;
end;

procedure TLjapunov.btnGrubLjaponvClick(Sender: TObject);
begin

end;

procedure TLjapunov.BtnLjapunov1Change(Sender: TObject);
Var
    i,j  : Integer;
  begin
  // A*P + P*A' = Q
   A_BKC (N,M,P,A,B,K,C,Ag);
   Az:=Ag; Qz:=V;
   lyap(N,Az,Qz,Py);
   btnQdif.Text:=FloatToStrF(QDif,ffNumber, 8, 4);
   StrGridLjapunov.ColCount:= N+1;
   StrGridLjapunov.RowCount:= N+1;
     //GridInit; //Beschriftung den Tabelle
     StrGridLjapunov.cells[0,0]:= 'Py';
     for i:=1 to N do StrGridLjapunov.cells[0,i]:=FloatToStr(i);
     for i:=1 to N do StrGridLjapunov.cells[i,0]:=FloatToStr(i);

   for j:=1 to N do begin
     for i:=1 to N do begin
        StrGridLjapunov.cells[j,i]:=FloatToStr(Py[j+N*(i-1)]);
     end;
   end;
end;

procedure TLjapunov.btnLjapunov2Change(Sender: TObject);
Var
    i,j  : Integer;
  begin
   // Ag'*P + P*Ag = Q
   A_BKC (N,M,P,A,B,K,C,Ag);
   QX_Aus (Q,C,K,R,Qx,N,M,P);
   Az:=Ag; Qz:=Qx;
   Transp(N,Az);

   lyap(N,Az,Qz,Px);
   btnQdif.Text:=FloatToStrF(QDif,ffNumber, 8, 4);
   StrGridLjapunov.ColCount:= N+1;
   StrGridLjapunov.RowCount:= N+1;
     //GridInit; //Beschriftung den Tabelle
     StrGridLjapunov.cells[0,0]:= 'Px';
     for i:=1 to N do StrGridLjapunov.cells[0,i]:=FloatToStr(i);
     for i:=1 to N do StrGridLjapunov.cells[i,0]:=FloatToStr(i);

   for j:=1 to N do begin
     for i:=1 to N do begin
        StrGridLjapunov.cells[j,i]:=FloatToStr(Px[j+N*(i-1)]);
     end;
   end;
  end;

procedure TLjapunov.btnMinvChange(Sender: TObject);
begin
  MaInvSel:=1;
end;

procedure TLjapunov.btnTmInvChange(Sender: TObject);
begin
  Mainvsel:=3;
end;

procedure TLjapunov.btnViewAChange(Sender: TObject);
Var
    i,j  : Integer;
begin
 StrGridLjapunov.ColCount:= N+1;
 StrGridLjapunov.RowCount:= N+1;
   //GridInit; //Beschriftung den Tabelle
   StrGridLjapunov.cells[0,0]:= 'A';
   for i:=1 to N do StrGridLjapunov.cells[0,i]:=FloatToStr(i);
   for i:=1 to N do StrGridLjapunov.cells[i,0]:=FloatToStr(i);

 for j:=1 to N do begin
   for i:=1 to N do begin
       StrGridLjapunov.cells[j,i]:=FloatToStr(A[j+N*(i-1)]);
   end;
 end;

end;

procedure TLjapunov.btnViewAtChange(Sender: TObject);
Var
    i,j  : Integer;
begin
 StrGridLjapunov.ColCount:= N+1;
 StrGridLjapunov.RowCount:= N+1;
   //GridInit; //Beschriftung den Tabelle
   StrGridLjapunov.cells[0,0]:= 'Ag';
   for i:=1 to N do StrGridLjapunov.cells[0,i]:=FloatToStr(i);
   for i:=1 to N do StrGridLjapunov.cells[i,0]:=FloatToStr(i);
  //  Mt:=A;
  // Transp(N,Mt);
 A_BKC (N,M,P,A,B,K,C,Ag);

 for j:=1 to N do begin
   for i:=1 to N do begin
       StrGridLjapunov.cells[j,i]:=FloatToStr(Ag[j+N*(i-1)]);
   end;
 end;

end;

procedure TLjapunov.btnViewEChange(Sender: TObject);
Var
    i,j  : Integer;
begin

 StrGridLjapunov.ColCount:= 3;   // Spalten
 StrGridLjapunov.RowCount:= N+1; // Zeilen
 StrGridLjapunov.cells[0,0]:= 'E Eig';
 StrGridLjapunov.cells[1,0]:= 'Real';
 StrGridLjapunov.cells[2,0]:= 'Imag';

 A_BKC (N,M,P,A,B,K,C,Ag);
 Ax:=Ag;
 EIGRF (N,Ax,U,Wr,Wi);

 for i:=1 to N do StrGridLjapunov.cells[0,i]:=FloatToStr(i);
 for i:=1 to N do StrGridLjapunov.cells[1,i]:=FloatToStr(Wr[i]);
 for i:=1 to N do StrGridLjapunov.cells[2,i]:=FloatToStr(Wi[i]);

end;

procedure TLjapunov.btnViewEvChange(Sender: TObject);
Var
    i,j  : Integer;
begin

 StrGridLjapunov.ColCount:= N+1;   // Spalten
 StrGridLjapunov.RowCount:= N+1; // Zeilen
 StrGridLjapunov.cells[0,0]:= 'E Vek';
 Ax:=A;
 EIGRF (N,Ax,TM,Wr,Wi);

 for i:=1 to N do StrGridLjapunov.cells[0,i]:=FloatToStr(i);
 for i:=1 to N do StrGridLjapunov.cells[i,0]:=FloatToStr(i);

  for j:=1 to N do begin
   for i:=1 to N do begin
      StrGridLjapunov.cells[j,i]:=FloatToStr(TM[j+N*(i-1)]);
   end;
 end;

end;

procedure TLjapunov.btnViewPChange(Sender: TObject);
Var
    i,j  : Integer;
begin
 StrGridLjapunov.ColCount:= N+1;
 StrGridLjapunov.RowCount:= N+1;
   //GridInit; //Beschriftung den Tabelle
   StrGridLjapunov.cells[0,0]:= 'P';
   for i:=1 to N do StrGridLjapunov.cells[0,i]:=FloatToStr(i);
   for i:=1 to N do StrGridLjapunov.cells[i,0]:=FloatToStr(i);

 for j:=1 to N do begin
   for i:=1 to N do begin
      StrGridLjapunov.cells[j,i]:=FloatToStr(Px[j+N*(i-1)]);
   end;
 end;


end;

procedure TLjapunov.btnViewQChange(Sender: TObject);
Var
    i,j  : Integer;
begin
 StrGridLjapunov.ColCount:= N+1;
 StrGridLjapunov.RowCount:= N+1;
   //GridInit; //Beschriftung den Tabelle
   StrGridLjapunov.cells[0,0]:= 'Q';
   for i:=1 to N do StrGridLjapunov.cells[0,i]:=FloatToStr(i);
   for i:=1 to N do StrGridLjapunov.cells[i,0]:=FloatToStr(i);

   QX_Aus (Q,C,K,R,Qx,N,M,P);

 for j:=1 to N do begin
   for i:=1 to N do begin
      StrGridLjapunov.cells[j,i]:=FloatToStr(Qx[j+N*(i-1)]);
   end;
 end;
end;

procedure TLjapunov.btnViewSUChange(Sender: TObject);
Var
    i,j  : Integer;
begin
 StrGridLjapunov.ColCount:= N+1;
 StrGridLjapunov.RowCount:= N+1;
   //GridInit; //Beschriftung den Tabelle
   StrGridLjapunov.cells[0,0]:= 'SU';
   for i:=1 to N do StrGridLjapunov.cells[0,i]:=FloatToStr(i);
   for i:=1 to N do StrGridLjapunov.cells[i,0]:=FloatToStr(i);

 for j:=1 to N do begin
   for i:=1 to N do begin
      StrGridLjapunov.cells[j,i]:=FloatToStr(SU[j+N*(i-1)]);
   end;
 end;

end;

procedure TLjapunov.btnViewUChange(Sender: TObject);
Var
    i,j  : Integer;
begin
 StrGridLjapunov.ColCount:= N+1;
 StrGridLjapunov.RowCount:= N+1;
   //GridInit; //Beschriftung den Tabelle
   StrGridLjapunov.cells[0,0]:= 'U';
   for i:=1 to N do StrGridLjapunov.cells[0,i]:=FloatToStr(i);
   for i:=1 to N do StrGridLjapunov.cells[i,0]:=FloatToStr(i);

 for j:=1 to N do begin
   for i:=1 to N do begin
      StrGridLjapunov.cells[j,i]:=FloatToStr(U[j+N*(i-1)]);
   end;
 end;

end;

procedure TLjapunov.btnQdifChange(Sender: TObject);
begin
  btnQdif.Text:=FloatToStrF(QDif,ffCurrency, 8, 4);
end;

procedure TLjapunov.FormCreate(Sender: TObject);
begin

end;

procedure TLjapunov.Label1Click(Sender: TObject);
begin

end;

procedure TLjapunov.btnJjapunovS2Change(Sender: TObject);
Var
    i,j  : Integer;
  begin
   // A'*P + P*A = Q
    A_BKC (N,M,P,A,B,K,C,Ag);
    QX_Aus (Q,C,K,R,Qx,N,M,P);
   Az:=Ag; Qz:=Qx;
   Transp(N,Az);
   LjpForm2:=False; Skip:=False;
   Ljapunow (N, Az,Qz,Pz,SU,LjpForm2,skip);
   N:=N;
   StrGridLjapunov.ColCount:= N+1;
   StrGridLjapunov.RowCount:= N+1;
     //GridInit; //Beschriftung den Tabelle
     StrGridLjapunov.cells[0,0]:= 'P';
     for i:=1 to N do StrGridLjapunov.cells[0,i]:=FloatToStr(i);
     for i:=1 to N do StrGridLjapunov.cells[i,0]:=FloatToStr(i);

   for j:=1 to N do begin
     for i:=1 to N do begin
        StrGridLjapunov.cells[j,i]:=FloatToStr(Pz[j+N*(i-1)]);
     end;
   end;

end;

procedure TLjapunov.btnJjapunovS1Change(Sender: TObject);
Var
    i,j  : Integer;
  begin
   // Ag'*P + P*Ag = V
   A_BKC (N,M,P,A,B,K,C,Ag);
   Az:=Ag; Qz:=V;

   LjpForm2:=False; Skip:=False;
   Ljapunow (N, Az,Qz,Pz,SU,LjpForm2,skip);
   N:=N;
   StrGridLjapunov.ColCount:= N+1;
   StrGridLjapunov.RowCount:= N+1;
     //GridInit; //Beschriftung den Tabelle
     StrGridLjapunov.cells[0,0]:= 'P';
     for i:=1 to N do StrGridLjapunov.cells[0,i]:=FloatToStr(i);
     for i:=1 to N do StrGridLjapunov.cells[i,0]:=FloatToStr(i);

   for j:=1 to N do begin
     for i:=1 to N do begin
        StrGridLjapunov.cells[j,i]:=FloatToStr(Pz[j+N*(i-1)]);
     end;
   end;
end;

procedure TLjapunov.btnMaInvChange(Sender: TObject);
begin
  Mainvsel:=2;
end;

procedure TLjapunov.rbBalanceChange(Sender: TObject);
Var
    i,j  : Integer;
  begin
   // Ag'*P + P*Ag = V
   A_BKC (N,M,P,A,B,K,C,Ag);
   Az:=Ag; Qz:=V;

   BALANCE (n,Ag,low, upp, Ew);

   StrGridLjapunov.ColCount:= N+1;
   StrGridLjapunov.RowCount:= N+1;
     //GridInit; //Beschriftung den Tabelle
     StrGridLjapunov.cells[0,0]:= 'Ag';
     for i:=1 to N do StrGridLjapunov.cells[0,i]:=FloatToStr(i);
     for i:=1 to N do StrGridLjapunov.cells[i,0]:=FloatToStr(i);

   for j:=1 to N do begin
     for i:=1 to N do begin
        StrGridLjapunov.cells[j,i]:=FloatToStr(Ag[j+N*(i-1)]);
     end;
   end;



end;

procedure TLjapunov.rbEwHesChange(Sender: TObject);
Var
    i,j  : Integer;
  begin
  
 StrGridLjapunov.ColCount:= 3;   // Spalten
 StrGridLjapunov.RowCount:= N+1; // Zeilen
 StrGridLjapunov.cells[0,0]:= 'EwHes';
 StrGridLjapunov.cells[1,0]:= 'Real';
 StrGridLjapunov.cells[2,0]:= 'Imag';

 for i:=1 to N do StrGridLjapunov.cells[0,i]:=FloatToStr(i);
 for i:=1 to N do StrGridLjapunov.cells[1,i]:=FloatToStr(Ew[i]);
// for i:=1 to N do StrGridLjapunov.cells[2,i]:=FloatToStr(Wi[i]);


end;

procedure TLjapunov.rbQrthesChange(Sender: TObject);
Var
    i,j  : Integer;
  begin

   A_BKC (N,M,P,A,B,K,C,Ag);
   Az:=Ag; Qz:=V;

   BALANCE (n,Ag,low, upp, Ew);

   Orthes (n,low, upp, Ag, Ew);

   StrGridLjapunov.ColCount:= N+1;
   StrGridLjapunov.RowCount:= N+1;
     //GridInit; //Beschriftung den Tabelle
     StrGridLjapunov.cells[0,0]:= 'Ag';
     for i:=1 to N do StrGridLjapunov.cells[0,i]:=FloatToStr(i);
     for i:=1 to N do StrGridLjapunov.cells[i,0]:=FloatToStr(i);

   for j:=1 to N do begin
     for i:=1 to N do begin
        StrGridLjapunov.cells[j,i]:=FloatToStr(Ag[j+N*(i-1)]);
     end;
   end;

end;

procedure TLjapunov.rbQrTransChange(Sender: TObject);
Var
    i,j  : Integer;
  begin

   A_BKC (N,M,P,A,B,K,C,Ag);
   Az:=Ag; Qz:=V;

   BALANCE (n,Ag,low, upp, Ew);

   Orthes (n,low, upp, Ag, Ew);

   ORTRANS (n,low,upp, Ag, He, Ew);

   StrGridLjapunov.ColCount:= N+1;
   StrGridLjapunov.RowCount:= N+1;
     //GridInit; //Beschriftung den Tabelle
     StrGridLjapunov.cells[0,0]:= 'Hess';
     for i:=1 to N do StrGridLjapunov.cells[0,i]:=FloatToStr(i);
     for i:=1 to N do StrGridLjapunov.cells[i,0]:=FloatToStr(i);

   for j:=1 to N do begin
     for i:=1 to N do begin
        StrGridLjapunov.cells[j,i]:=FloatToStr(He[j+N*(i-1)]);
     end;
   end;

end;

end.

