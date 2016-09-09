{..............................................................................}
{                   }
{                                        }
{..............................................................................}

{..............................................................................}

//Global Constants
Const
	InputFileName = 'C:\Users\admin\Documents\GitHub\PcbTools\PCBTools\AutoPowerPlaneCut\Vertices.txt';
	TrackWidth = 15.0;
	PlaneCutLayer = eInternalPlane2;
	
Var
	Board     : IPCB_Board;
    ItemCount : Integer;
	
Procedure CreateATrack(X1, Y1, X2, Y2, Width, Net, LayerID);
Var
    Track       : IPCB_Track;
Begin
    Track          := PCBServer.PCBObjectFactory(eTrackObject, eNoDimension, eCreate_Default);
    PCBServer.SendMessageToRobots(Track.I_ObjectAddress ,c_Broadcast, PCBM_BeginModify ,c_NoEventData);
    Track.Width    := MilsToCoord(Width);
    Track.Layer    := LayerID;
    //Track.Net    := Net;
    Track.X1       := MilsToCoord(X1) + Board.XOrigin;
    Track.Y1       := MilsToCoord(Y1) + Board.YOrigin;
    Track.X2       := MilsToCoord(X2) + Board.XOrigin;
    Track.Y2       := MilsToCoord(Y2) + Board.YOrigin;
    Board.AddPCBObject(Track); 
    Track.Selected := True;
    PCBServer.SendMessageToRobots(Track.I_ObjectAddress, c_Broadcast, PCBM_EndModify, c_NoEventData);
    PCBServer.SendMessageToRobots(Board.I_ObjectAddress, c_Broadcast, PCBM_BoardRegisteration, Track.I_ObjectAddress);
End;

//extract the name of the board from the board fileName (full path)
Function Split(MyString, Delimeter) : TStrings;
var
    List    : TStrings;
    aString : Tstring;
Begin
    List               := TStringList.Create;
    //sepearate the board's name according to file path slashes
    List.Delimiter     := Delimeter;  
    List.DelimitedText := MyString;
    //last item in the list will be FileName.PCBDOC
    Result             := List;
End;

Function Clamp(Input, Lower, Upper) : Float;
Begin
	If Input < Lower Then Result := Lower
	Else if Input > Upper Then Result := Upper
	Else Result := Input;
		
End;
		
Procedure AddCutsFromFile(Dummy);
var
    InputFile  : TextFile;
    I          : Integer;
    Line       : String;
	X1, Y1       : Float;
	X2, Y2       : Float;
	LineSeparated : TStrings;
	
	IndP0 : Integer;
	IndP1 : Integer;
	
Begin
   
    AssignFile(InputFile, InputFileName);
    Reset(InputFile);

    Try
        while not EOF(InputFile) do
        begin
            readln(InputFile, Line);
            If Not VarIsNull(Line) Then
            Begin
				LineSeparated := Split(Line, ' ');
				For I := 0 to (LineSeparated.Count/2 - 1) Do
				Begin
					IndP0 := 2*I;
					IndP1 := 2*(I+1) mod LineSeparated.Count;
					X1 := StrToFloat(LineSeparated[IndP0]);
					Y1 := StrToFloat(LineSeparated[IndP0+1]);
					X2 := StrToFloat(LineSeparated[IndP1]);
					Y2 := StrToFloat(LineSeparated[IndP1+1]);
					CreateATrack(X1, Y1, X2, Y2, TrackWidth, 0, PlaneCutLayer);
				End;
				LineSeparated.Free;
				Inc(ItemCount);
				//Break;
				
            End
        end;

    Finally
        CloseFile(InputFile);
    End;
End;

Procedure ProcessPCB;
Var

    Item      : IPCB_Primitive;
    Iterator  : IPCB_BoardIterator;
   
    NetName   : String;
    HoleSize  : String;
    PadSize   : String;
    CenterX   : String;
    CenterY   : String;
    ItemsList : TStringList;
    ItemString : String;
    NetsToExport : TStringList;

Begin
    ItemCount       := 0;

    // Retrieve the current board
    Board := PCBServer.GetCurrentPCBBoard;
    If Board = Nil Then Exit;
	
	AddCutsFromFile(0);
    
	

    // Display the count result on a dialog.
    ShowMessage('Item Count = ' + IntToStr(ItemCount));
End;
{..............................................................................}

{..............................................................................}
