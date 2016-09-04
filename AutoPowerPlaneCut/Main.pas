{..............................................................................}
{                   }
{                                        }
{..............................................................................}

{..............................................................................}

//var
//   Board      : IPCB_Board;
   

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

Begin
    ItemCount       := 0;

    // Retrieve the current board
    Board := PCBServer.GetCurrentPCBBoard;
    If Board = Nil Then Exit;
    
    ItemsList := TStringList.Create;

    // retrieve the iterator
    Iterator        := Board.BoardIterator_Create;
    Iterator.AddFilter_ObjectSet(MkSet(ePadObject, eViaObject));
    Iterator.AddFilter_LayerSet(MkSet(eMultiLayer));
    Iterator.AddFilter_Method(eProcessAll);

    // Search and count pads
    Item := Iterator.FirstPCBObject;
    While (Item <> Nil) Do
    Begin
        Inc(ItemCount);
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
			ShowMessage(ItemString);
			
			ItemsList.Add(ItemString);
		End;
    
        Item := Iterator.NextPCBObject;
    End;
    Board.BoardIterator_Destroy(Iterator);
    WriteOutput('ItemsList.txt', ItemsList);
    ItemsList.Free;

    // Display the count result on a dialog.
    ShowMessage('Item Count = ' + IntToStr(ItemCount));
End;
{..............................................................................}

{..............................................................................}
