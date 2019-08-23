

--#########################################################
function TextBoxClass.ProcessEvents(this)
    if (WM_MOUSEMOVE.occured) then
        if     (WM_MOUSEMOVE.iPar5==this.obj1Num) then
            chiWindowSetCursor(3);
        end
    end
    
    if (WM_LBUTTONDOWN.occured) then
        if     (WM_LBUTTONDOWN.iPar5==this.obj1Num) then
            this.OnSelection(this);
        else
            this.OnDeselection(this);
        end
    end
    
    if (WM_KEYDN.occured) then
        this.KeyDown(this);
    end
    
    if (WM_SHIFTBUTTONDN.occured) then
        this.shiftDown=true;
    end
    
    if (WM_SHIFTBUTTONUP.occured) then
        this.shiftDown=false;
    end
    
    if (WM_CHAR.occured) then
        this.KeyPress(this);
        if (not (this.CustomKeyPress==nil)) then
            this.CustomKeyPress(this);
        end
    end
    
    this.BlinkCursor(this)
    
    for k=1,this.eventCallbackCount do
        this.eventCallbacks[k](this);
    end
end

--#########################################################
function TextBoxClass.SetProperty(this,property,value)
    if     (property=="Master") then
        this.master=value;
        
        value.slaveCount=value.slaveCount+1;
        k=value.slaveCount
        value.slaves[k]=this;
    elseif (property=="Float") then
        this.float=value;
    elseif (property=="Text") then
        this.text=value;
        r=this.fontColor[1];
        g=this.fontColor[2];
        b=this.fontColor[3];
        a=this.fontColor[4];
        chiSetLabel3D(this.textNum,this.text,this.xpos+this.marginLeft+this.paddingLeft,this.ypos+this.marginBot+this.paddingBot,r,g,b,a,this.fontType)
    elseif (property=="Color") then
        this.fontColor[1]=value[1];
        this.fontColor[2]=value[2];
        this.fontColor[3]=value[3];
        this.fontColor[4]=value[4];
        r=value[1];
        g=value[2];
        b=value[3];
        a=value[4];
        chiSetLabel3D(this.textNum,this.text,this.xpos+this.marginLeft+this.paddingLeft,this.ypos+this.marginBot+this.paddingBot,r,g,b,a,this.fontType)
     elseif (property=="ShowOutline") then  
        if (value) then
            chi3DLineChangeColor(this.lineNum,col[1],col[2],col[3],col[4]);
        else
            chi3DLineChangeColor(this.lineNum,col[1],col[2],col[3],0.0);
        end
    end
end

--#########################################################
function TextBoxClass.Redraw(this)
    chiTransformSetScale(this.transform1,this.xSize-1,this.ySize-1,1.0)
    chiTransformSetTranslation(this.transform1,this.xpos+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth)

    chi3DLineChangeVertex(this.lineNum,0,this.xpos+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth);
    chi3DLineChangeVertex(this.lineNum,1,this.xpos+this.paddingLeft,this.ypos+this.ySize+this.paddingBot,this.zDepth);
    chi3DLineChangeVertex(this.lineNum,2,this.xpos+this.xSize+this.paddingLeft,this.ypos+this.ySize+this.paddingBot,this.zDepth);
    chi3DLineChangeVertex(this.lineNum,3,this.xpos+this.xSize+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth);
    chi3DLineChangeVertex(this.lineNum,4,this.xpos-1+this.paddingLeft,this.ypos+this.paddingBot,this.zDepth);
    
    if (this.inputBox) then
    
        if (this.selected) then
            r=this.fontColor[1];
            g=this.fontColor[2];
            b=this.fontColor[3];
            a=this.fontColor[4];
            chiSetLabel3D(this.textNum,this.text,this.xpos+this.marginLeft+this.paddingLeft,this.ypos+this.marginBot+this.paddingBot,r,g,b,a,this.fontType)
            chiSetLabel3D(this.inputNum,this.inputText,this.xpos+this.marginLeft+this.paddingLeft,this.ypos+this.marginBot+this.paddingBot,r,g,b,0,this.fontType)
            
        else
        
            if (this.text=="") then
                r=this.inputfontColor[1];
                g=this.inputfontColor[2];
                b=this.inputfontColor[3];
                a=this.inputfontColor[4];
                chiSetLabel3D(this.inputNum,this.inputText,this.xpos+this.marginLeft+this.paddingLeft,this.ypos+this.marginBot+this.paddingBot,r,g,b,a,this.fontType)
            
            else
                r=this.fontColor[1];
                g=this.fontColor[2];
                b=this.fontColor[3];
                a=this.fontColor[4];
                chiSetLabel3D(this.textNum,this.text,this.xpos+this.marginLeft+this.paddingLeft,this.ypos+this.marginBot+this.paddingBot,r,g,b,a,this.fontType)
                chiSetLabel3D(this.inputNum,this.inputText,this.xpos+this.marginLeft+this.paddingLeft,this.ypos+this.marginBot+this.paddingBot,r,g,b,0,this.fontType)
                --print(this.name,"MUFF",this.xpos+this.marginLeft+this.paddingLeft)
            end
        end
        chiSetLabel3D(this.selectNum,this.text,this.xpos+this.marginLeft+this.paddingLeft,this.ypos+this.marginBot+this.paddingBot,r,g,b,0,this.fontType)
    else
    
        r=this.fontColor[1];
        g=this.fontColor[2];
        b=this.fontColor[3];
        a=this.fontColor[4];
        chiSetLabel3D(this.textNum,this.text,this.xpos+this.marginLeft+this.paddingLeft,this.ypos+this.marginBot+this.paddingBot,r,g,b,a,this.fontType)
        
        if (this.selected) then
            r=this.selectfontColor[1];
            g=this.selectfontColor[2];
            b=this.selectfontColor[3];
            a=this.selectfontColor[4];
            chiSetLabel3D(this.selectNum,this.text,this.xpos+this.marginLeft+this.paddingLeft,this.ypos+this.marginBot+this.paddingBot,r,g,b,a,this.fontType)
        else
            chiSetLabel3D(this.selectNum,this.text,this.xpos+this.marginLeft+this.paddingLeft,this.ypos+this.marginBot+this.paddingBot,r,g,b,0,this.fontType)
        end
        chiSetLabel3D(this.inputNum,this.inputText,this.xpos+this.marginLeft+this.paddingLeft,this.ypos+this.marginBot+this.paddingBot,r,g,b,0,this.fontType)
    end   
    
    chiLineChangeVertex(this.cursorlineNum,0,this.xpos+this.paddingLeft+this.marginLeft+this.cursorX,this.ypos+this.paddingBot+2,0.0);
    chiLineChangeVertex(this.cursorlineNum,1,this.xpos+this.paddingLeft+this.marginLeft+this.cursorX,this.ypos+this.paddingBot+this.ySize-this.paddingTop-3,0.0);
    
end




--#########################################################
function TextBoxClass.SizeChanged(this)
    if (this.master~=nil) then
        --print(this.master.cursorY)
        if (not this.float) then
            this.xpos=this.master.cursorX;
            this.ypos=this.master.cursorY-this.ySize;
           
            chiSetLabelProperty3D(this.textNum,"Viewport",this.master.xmin,this.master.ymin,this.master.xmax,this.master.ymax);
            chiViewportSetProperty(this.obj1View,this.master.xmin+1,this.master.ymin+1,this.master.xmax,this.master.ymax);
            this.Redraw(this)
        end
        
        

    end
end

--#########################################################
function TextBoxClass.OnSelection(this)
    this.selected=true;
    --Compute cursor position
    textLen=string.len(this.text);
    curCursorPos=this.xpos+this.paddingLeft+this.marginLeft;
    preCursorPos=this.xpos+this.paddingLeft+this.marginLeft;
    this.cursorPosition=0;
    if (textLen==0) then
        finCursorPos=curCursorPos;
    else
        for k=1,textLen do
            keycode=string.byte(string.sub(this.text,k,k));
            curCursorPos=curCursorPos+chiGetCharacterWidth(keycode,this.fontType);
            if ((WM_LBUTTONDOWN.iPar0>=preCursorPos) and (WM_LBUTTONDOWN.iPar0<curCursorPos)) then
                curCursorPos = preCursorPos;
                this.curCursorPos=k;
                break;
            end
            this.cursorPosition=k;
            preCursorPos=preCursorPos+chiGetCharacterWidth(keycode,this.fontType);
        end
        finCursorPos=curCursorPos;
    end
    this.cursorX=finCursorPos-(this.xpos+this.paddingLeft+this.marginLeft);
    this.Redraw(this);
end

--#########################################################
function TextBoxClass.OnDeselection(this)
    this.selected=false;
    this.Redraw(this);
end

--#########################################################
function TextBoxClass.BlinkCursor(this)
    if (this.selected) then
        if (chi_programTime>(this.lastBlinkUpdated+1.0)) then
            this.lastBlinkUpdated=chi_programTime;
            if (this.cursorShown) then
                this.cursorShown=false;
                chiLineChangeColor(this.cursorlineNum,0.0,0.0,0.0,0.0);
            else
                this.cursorShown=true;
                chiLineChangeColor(this.cursorlineNum,0.0,0.0,0.0,1.0);
            end
        end
    else
        chiLineChangeColor(this.cursorlineNum,0.0,0.0,0.0,0.0);
    end
    
end

--#########################################################
function TextBoxClass.KeyDown(this)

    if (this.selected) then
    
        if (WM_KEYDN.iPar3==39) then
            
            keycode=string.byte(this.text,this.cursorPosition+1)
            displacement=chiGetCharacterWidth(keycode,this.fontType);
            if (displacement>0) then
                this.cursorPosition=this.cursorPosition+1;
                this.cursorX=this.cursorX+displacement;
                
                this.Redraw(this)
            end
        end
        if (WM_KEYDN.iPar3==37) then
            if (this.cursorPosition>0) then
                keycode=string.byte(this.text,this.cursorPosition)
                displacement=chiGetCharacterWidth(keycode,this.fontType);
            else
                displacement=0;
            end
            if (displacement>0) then
                this.cursorPosition=this.cursorPosition-1;
                this.cursorX=this.cursorX-displacement;
                
                this.Redraw(this)
            end
        end
        if (WM_KEYDN.iPar3==46) then
            temp=string.sub(this.text,1,this.cursorPosition);
            temp=temp..string.sub(this.text,this.cursorPosition+2);
            this.text=temp;

            this.Redraw(this)
        end
    end
end

--#########################################################
function TextBoxClass.KeyPress(this)
    if (this.selected) then
        keycode=WM_CHAR.iPar0;
        if (keycode~=8) then

            displacement=chiGetCharacterWidth(keycode,this.fontType);
            if (displacement>0) then
                
                
                temp=string.sub(this.text,1,this.cursorPosition);
                temp=temp..string.char(keycode)
                temp=temp..string.sub(this.text,this.cursorPosition+1);
                this.text=temp;
                this.cursorPosition=this.cursorPosition+1;
                this.cursorX=this.cursorX+displacement;
                
                this.Redraw(this)
            end
        else
            keycode=string.byte(this.text,this.cursorPosition)
            if (this.cursorPosition>1) then
                temp=string.sub(this.text,1,this.cursorPosition-1);
            else
                temp="";
            end
            temp=temp..string.sub(this.text,this.cursorPosition+1);
            this.text=temp;
            if ((this.cursorPosition>0) and (this.cursorX<this.xSize)) then
                this.cursorPosition=this.cursorPosition-1;
                this.cursorX=this.cursorX-chiGetCharacterWidth(keycode,this.fontType);
            end

            this.Redraw(this)
        end
    end
end

--#########################################################
function TextBoxClass.Hide(this)
    this.selected=false;
    this.SetProperty(this,"Color",{0,0,0,0});
    chiObjectSetProperty(this.obj1Num,"Hidden",true);
    this.hidden=true;
end

--#########################################################
function TextBoxClass.UnHide(this)
    r=0;
    g=0;
    b=0;
    a=1;
    this.SetProperty(this,"Color",{r,g,b,a});
    chiObjectSetProperty(this.obj1Num,"Hidden",false);
    this.hidden=false;
end