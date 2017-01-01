{..............................................................................}
{                   }
{                                        }
{..............................................................................}

{..............................................................................}

//var
//   Board      : IPCB_Board;
   

Function GetNetsForExport() : TStringList;
Var
    List : TStringList;
Begin
    List := TStringList.Create;
    
    //List.Add('SGND');
    //List.Add('VOUT_REG');
    //List.Add('VOUT_SENSOR');
    //List.Add('VDD_3.3');
    List.Add('+1V8');
    List.Add('+1V2');
    List.Add('VDD_MCU');
    List.Add('VREF+');
    List.Add('VBUS_HS');
    List.Add('VBUS_FS');
    List.Add('U5V_ST_LINK');
    List.Add('EMU_5V');
    List.Add('EMU_3V3');
    List.Add('E5V');
    Result := List
End;

Function NetInList(NetName, NetsToExport) : Boolean;
Var
    i            : Integer;

Begin
    //NetsToExport := GetNetsForExport();
    Result := False;
    
    For i := 0 to NetsToExport.Count - 1 do
        If NetName = NetsToExport[i] Then
        Begin
            Result := True;
            Exit;
        End;
    
End;
   
   
Function Max(A, B) : Int;
Begin
    If A > B Then
    Begin
        Result := A;
    End
    Else
    Begin
        Result := B;
    End;
End;



Procedure WriteOutput(FileName, ItemsList);
Var
    OutputPath : String;

Begin
    
    OutputPath := ExtractFilePath(PCBServer.GetCurrentPCBBoard.FileName);

    ItemsList.SaveToFile(OutputPath + FileName)

End;

Procedure AddBoardOutlineToList(PCB_Board; ItemsList);
Var
    I,J       : Integer;
    X,Y       : Float;
    ItemString : String;
Begin
    PCB_Board.BoardOutline.Invalidate;
    PCB_Board.BoardOutline.Rebuild;
    PCB_Board.BoardOutline.Validate;

    ItemsList.Add('***BOARD OUTLINE***');
    //PCBServer.PreProcess;
    // Step through each of the vertices of the Board Outline in turn.
    For I := 0 To PCB_Board.BoardOutline.PointCount - 1 Do
    Begin
        // Set the value of J to point to the "next" vertex; this is normally
        // I + 1, but needs to be set to 0 instead for the very last vertex
        // that is processed by this loop.
        If I = PCB_Board.BoardOutline.PointCount - 1 Then
            J := 0
        Else
            J := I + 1;

        If PCB_Board.BoardOutline.Segments[I].Kind = ePolySegmentLine Then
        Begin
            // Current segment is a straight line; create a Track object.

            X := FloatToStr(CoordToMils(PCB_Board.BoardOutline.Segments[I].vx - PCB_Board.XOrigin));
            Y := FloatToStr(CoordToMils(PCB_Board.BoardOutline.Segments[I].vy - PCB_Board.YOrigin));
            ItemString := 'TRACK' + ';' + X + ';' + Y + ';';
            X := FloatToStr(CoordToMils(PCB_Board.BoardOutline.Segments[J].vx - PCB_Board.XOrigin));
            Y := FloatToStr(CoordToMils(PCB_Board.BoardOutline.Segments[J].vy - PCB_Board.YOrigin));
            ItemString :=  ItemString + X + ';' + Y;
            ItemsList.Add(ItemString);
        End
        Else
        Begin
            // Current segment is an arc; create an Arc object.
            X := FloatToStr(CoordToMils(PCB_Board.BoardOutline.Segments[I].cx - PCB_Board.XOrigin));
            Y := FloatToStr(CoordToMils(PCB_Board.BoardOutline.Segments[I].cy - PCB_Board.YOrigin));
            ItemString := 'ARC' + ';' + X + ';' + Y + ';';

            ItemString := ItemString + FloatToStr(CoordToMils(PCB_Board.BoardOutline.Segments[I].Radius)) + ';';
            ItemString := ItemString + FloatToStr(PCB_Board.BoardOutline.Segments[I].Angle1) + ';';
            ItemString := ItemString + FloatToStr(PCB_Board.BoardOutline.Segments[I].Angle2);
        End;
    End;
    //PCBServer.PostProcess;

    // Display (unconditionally) the layer selected by the user.
    //PCB_Board.LayerIsDisplayed[ALayer] := True;

    // Refresh PCB workspace.
    //ResetParameters;
    //AddStringParameter('Action', 'Redraw');
    //RunProcess('PCB:Zoom');
End;

Procedure ProcessPCB;
Var
    Board     : IPCB_Board;
    Item      : IPCB_Primitive;
    Iterator  : IPCB_BoardIterator;
    ItemCount : Integer;
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
    
    ItemsList := TStringList.Create;
    NetsToExport := GetNetsForExport();

    // retrieve the iterator
    Iterator        := Board.BoardIterator_Create;
    Iterator.AddFilter_ObjectSet(MkSet(ePadObject, eViaObject));
    Iterator.AddFilter_LayerSet(MkSet(eMultiLayer));
    Iterator.AddFilter_Method(eProcessAll);

    // Search and count pads
    Item := Iterator.FirstPCBObject;
    While (Item <> Nil) Do
    Begin
        {
        If Item.Net <> Nil Then
        Begin
            NetName := Item.Net.Name;
        End
        Else
        Begin
            NetName := 'No Net';
        End;
        }
        
        If Item.Net <> Nil Then
        Begin
            NetName := Item.Net.Name;
            If NetInList(NetName, NetsToExport) Then
            Begin
                Inc(ItemCount);
            
                HoleSize := FloatToStr(CoordToMils(Item.HoleSize));
            
                If Item.ObjectId = eViaObject Then
                Begin
                    PadSize := FloatToStr(CoordToMils(Item.Size));
                End
                Else If Item.ObjectId = ePadObject Then
                Begin
                    PadSize := FloatToStr(CoordToMils(Max(Item.MidXSize, Item.MidYSize)));
                End;

                CenterX := FloatToStr(CoordToMils(Item.X - Board.XOrigin));
                CenterY := FloatToStr(CoordToMils(Item.Y - Board.YOrigin));

                ItemString := NetName + ';' + HoleSize + ';' + PadSize + ';' + CenterX + ';' + CenterY;
                //ShowMessage(ItemString);

                ItemsList.Add(ItemString);
            End;
        End;
    
        Item := Iterator.NextPCBObject;
    End;
    Board.BoardIterator_Destroy(Iterator);


    AddBoardOutlineToList(Board, ItemsList);

    WriteOutput('ItemsList.txt', ItemsList);
    ItemsList.Free;

    // Display the count result on a dialog.
    ShowMessage('Item Count = ' + IntToStr(ItemCount));
End;
{..............................................................................}

{..............................................................................}
