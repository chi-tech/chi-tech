--=============================================== Console


consoleInputBox=TextBoxClass.New("ConsoleInputBox"..99);

consoleInputBox.inputBox=true;
consoleInputBox.paddingTop=0
consoleInputBox.paddingLeft=5;
consoleInputBox.paddingBot=-5;
consoleInputBox.xSize=600;
consoleInputBox.ySize=20;
consoleInputBox.text="";
consoleInputBox.inputText="Type lua command";
consoleInputBox.SetProperty(consoleInputBox,"Master",panels[6])
consoleInputBox.selected=true

consoleButton=ButtonClass.New("Run")
consoleButton.SetProperty(consoleButton,"Master",panels[6])
consoleButton.float=false;
consoleButton.paddingBot=0;
consoleButton.ySize=18;
consoleButton.yoffset=-7;
consoleButton.textPadding=4;
dRegisterFeature(consoleButton);

--Execute console command
function consoleInputBox.CustomKeyPress(this)
    if (WM_CHAR.iPar0==13) then
        chunk,err=load(this.text)
        if (chunk==nil) then
            print(err);
        else
            chunk()
        end
        this.text="";
        this.cursorPosition=0;
        this.cursorX=0;
        this.Redraw(this);
    end
end

consoleInputBox.SetProperty(consoleInputBox,"Float",false);
--consoleInputBox.text="propertyGrid.SetProperty(propertyGrid,\"Master\",panels[5])"
dRegisterFeature(consoleInputBox);

console=ConsoleClass.New("Main Console");
console.SetProperty(console,"Master",panels[6])
dRegisterFeature(console);