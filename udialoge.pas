unit uDialoge;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, Menus, ComCtrls,
  StdCtrls, Grids, ExtCtrls, LazUnicode, uTools, uLak, uDaten, uViewInput, urunmode;

type

  { TOptArf }

  TOptArf = class(TForm)
    MainMenu1: TMainMenu;
    MenuDatei: TMenuItem;
    MenuIInpView: TMenuItem;
    MenuLjapunov: TMenuItem;
    MenuOptMode: TMenuItem;
    MenuOptionen: TMenuItem;
    MenuIExe: TMenuItem;
    MenuLoad: TMenuItem;
    MenuExit: TMenuItem;
    OpenDialog1: TOpenDialog;
    procedure FormCreate(Sender: TObject);
    procedure MenuExitClick(Sender: TObject);
    procedure MenuIInpViewClick(Sender: TObject);
    procedure MenuLjapunovClick(Sender: TObject);
    procedure MenuLoadClick(Sender: TObject);
    procedure MenuOptionenClick(Sender: TObject);
  private

  public

    Var
     OptArf : TOptArf;

  end;



implementation

{$R *.lfm}

{ TOptArf }

procedure TOptArf.MenuLoadClick(Sender: TObject);

var Filename : String;
    DimMat   : String;
    RealDim  : real;

begin
 if OpenDialog1.Execute then
  begin
    RealDim :=SizeOf(RealDim);
        fmViewInput.show;
    Filename := OpenDialog1.FileName;
  //  ShowMessage(Filename);

    loadfile ( FileName, A, B, C, K, Q, R, V, N, M, P); // M-File

    fmViewInput.MatDimN.Text:= IntToStr(N);
    fmViewInput.MatDimM.Text:= IntToStr(M);
    fmViewInput.MatDimP.Text:= IntToStr(P);

  end;


end;

procedure TOptArf.MenuOptionenClick(Sender: TObject);
begin

end;

procedure TOptArf.MenuExitClick(Sender: TObject);
begin
  Close;
end;

procedure TOptArf.MenuIInpViewClick(Sender: TObject);
begin
  fmViewInput.show;
end;

procedure TOptArf.MenuLjapunovClick(Sender: TObject);
begin
  Ljapunov.show;
end;

procedure TOptArf.FormCreate(Sender: TObject);
begin
 N:=1;M:=1;P:=1;
 A[1]:=-1;B[1]:=1;C[1]:=1;K[1]:=1;Q[1]:=1;R[1]:=1;V[1]:=1;

end;

end.

