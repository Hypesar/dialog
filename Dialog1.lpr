program Dialog1;

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}
  cthreads,
  {$ENDIF}
  {$IFDEF HASAMIGA}
  athreads,
  {$ENDIF}
  Interfaces, // this includes the LCL widgetset
  Forms, uDialoge, uViewInput, uDaten, urunmode, utools
  { you can add units after this };

{$R *.res}

var OptArf : TOptArf;

begin
  RequireDerivedFormResource:=True;
  Application.Scaled:=True;
  Application.Initialize;
  Application.CreateForm(TOptArf, OptArf);
  Application.CreateForm(TfmViewInput, fmViewInput);
  Application.CreateForm(TLjapunov, Ljapunov);
  Application.Run;
end.


