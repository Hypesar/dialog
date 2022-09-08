unit uViewInput;

{$mode ObjFPC}{$H+}

interface

uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, Menus, StdCtrls, Grids,
  uDaten, ExtCtrls;

type

  { TfmViewInput }

  TfmViewInput = class(TForm)
    MatDimN: TEdit;
    MatDimM: TEdit;
    MatDimP: TEdit;
    MainMenu1: TMainMenu;
    MenuItem1: TMenuItem;
    MenuItem2: TMenuItem;
    btnA: TRadioButton;
    btnB: TRadioButton;
    BtnC: TRadioButton;
    btnK: TRadioButton;
    btnV: TRadioButton;
    btnQ: TRadioButton;
    btnR: TRadioButton;
    StaticText1: TStaticText;
    StringGrid1: TStringGrid;
    procedure btnBChange(Sender: TObject);
    procedure BtnCChange(Sender: TObject);
    procedure btnVChange(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure MatDimNChange(Sender: TObject);
    procedure MenuItem1Click(Sender: TObject);
    procedure MenuItem2Click(Sender: TObject);
    procedure btnAChange(Sender: TObject);
    procedure btnKChange(Sender: TObject);
    procedure btnQChange(Sender: TObject);
    procedure btnRChange(Sender: TObject);
  private

  public

  end;

var
  fmViewInput: TfmViewInput;

implementation

{$R *.lfm}

{ TfmViewInput }

procedure TfmViewInput.FormCreate(Sender: TObject);
begin
 fmViewInput.MatDimN.Text:= IntToStr(N);
 fmViewInput.MatDimM.Text:= IntToStr(M);
 fmViewInput.MatDimP.Text:= IntToStr(P);

end;

procedure TfmViewInput.btnBChange(Sender: TObject);
var
  i,j : Integer;
begin
  StringGrid1.caption := 'B';
  StringGrid1.ColCount:= M+1;   // Spalten
  StringGrid1.RowCount:= N+1;   // Zeilen
  //GridInit; //Beschriftung den Tabelle
   fmViewInput.StringGrid1.cells[0,0]:= 'B';
   for i:=1 to N do fmViewInput.StringGrid1.cells[0,i]:=FloatToStr(i);
   for i:=1 to M do fmViewInput.StringGrid1.cells[i,0]:=FloatToStr(i);

 for j:=1 to N do begin
   for i:=1 to M do begin
      fmViewInput.StringGrid1.cells[i,j]:=FloatToStr(B[i+M*(j-1)]);
   end;
 end;
end;

procedure TfmViewInput.BtnCChange(Sender: TObject);
var
  i,j  : Integer;
begin
 StringGrid1.ColCount:= N+1;
 StringGrid1.RowCount:= P+1;
   //GridInit; //Beschriftung den Tabelle
   fmViewInput.StringGrid1.cells[0,0]:= 'C';
   for i:=1 to P do fmViewInput.StringGrid1.cells[0,i]:=FloatToStr(i);
   for i:=1 to N do fmViewInput.StringGrid1.cells[i,0]:=FloatToStr(i);

 for j:=1 to P do begin
  for i:=1 to N do begin
   fmViewInput.StringGrid1.cells[i,j]:=FloatToStr(C[i+N*(j-1)]);
  end;
 end;

end;

procedure TfmViewInput.btnVChange(Sender: TObject);

var
  i,j  : Integer;
begin
 StringGrid1.ColCount:= N+1;
 StringGrid1.RowCount:= N+1;
   //GridInit; //Beschriftung den Tabelle
   fmViewInput.StringGrid1.cells[0,0]:= 'V';
   for i:=1 to N do fmViewInput.StringGrid1.cells[0,i]:=FloatToStr(i);
   for i:=1 to N do fmViewInput.StringGrid1.cells[i,0]:=FloatToStr(i);

 for j:=1 to N do begin
  for i:=1 to N do begin
   fmViewInput.StringGrid1.cells[j,i]:=FloatToStr(V[j+N*(i-1)]);
  end;
 end;
end;

procedure TfmViewInput.MatDimNChange(Sender: TObject);

begin

end;

procedure TfmViewInput.MenuItem1Click(Sender: TObject);
begin

end;

procedure TfmViewInput.MenuItem2Click(Sender: TObject);
begin
  close;
end;

procedure TfmViewInput.btnAChange(Sender: TObject);

var
  i,j  : Integer;
begin
 StringGrid1.ColCount:= N+1;
 StringGrid1.RowCount:= N+1;
 //GridInit; //Beschriftung den Tabelle
  fmViewInput.StringGrid1.cells[0,0]:= 'A';
  for i:=1 to N do fmViewInput.StringGrid1.cells[0,i]:=FloatToStr(i);
  for i:=1 to N do fmViewInput.StringGrid1.cells[i,0]:=FloatToStr(i);

 for j:=1 to N do begin
  for i:=1 to N do begin
   fmViewInput.StringGrid1.cells[j,i]:=FloatToStr(A[j+N*(i-1)]);
  end;
 end;
end;

procedure TfmViewInput.btnKChange(Sender: TObject);

var
  i,j  : Integer;
begin
  StringGrid1.ColCount:= M+1;
  StringGrid1.RowCount:= P+1;
  //GridInit; //Beschriftung den Tabelle
   fmViewInput.StringGrid1.cells[0,0]:= 'K';
   for i:=1 to M do fmViewInput.StringGrid1.cells[0,i]:=FloatToStr(i);
   for i:=1 to P do fmViewInput.StringGrid1.cells[i,0]:=FloatToStr(i);

    for j:=1 to M do begin
     for i:=1 to P do begin

      fmViewInput.StringGrid1.cells[j,i]:=FloatToStr(K[j+P*(i-1)]);

   //GridInit; //Beschriftung den Tabelle
     end;
    end;
end;

procedure TfmViewInput.btnQChange(Sender: TObject);

var
  i,j  : Integer;
begin
 StringGrid1.ColCount:= N+1;
 StringGrid1.RowCount:= N+1;
 //GridInit; //Beschriftung den Tabelle
  fmViewInput.StringGrid1.cells[0,0]:= 'Q';
  for i:=1 to N do fmViewInput.StringGrid1.cells[0,i]:=FloatToStr(i);
  for i:=1 to N do fmViewInput.StringGrid1.cells[i,0]:=FloatToStr(i);

 for j:=1 to N do begin
  for i:=1 to N do begin

   fmViewInput.StringGrid1.cells[j,i]:=FloatToStr(Q[j+N*(i-1)]);

//GridInit; //Beschriftung den Tabelle
  end;
 end;
end;

procedure TfmViewInput.btnRChange(Sender: TObject);

var
  i,j  : Integer;
begin
 StringGrid1.ColCount:= M+1;
 StringGrid1.RowCount:= M+1;
 //GridInit; //Beschriftung den Tabelle
  fmViewInput.StringGrid1.cells[0,0]:= 'R';
  for i:=1 to M do fmViewInput.StringGrid1.cells[0,i]:=FloatToStr(i);
  for i:=1 to M do fmViewInput.StringGrid1.cells[i,0]:=FloatToStr(i);

 for j:=1 to M do begin
  for i:=1 to M do begin

   fmViewInput.StringGrid1.cells[j,i]:=FloatToStr(R[j+M*(i-1)]);

//GridInit; //Beschriftung den Tabelle
  end;
 end;
end;

end.

