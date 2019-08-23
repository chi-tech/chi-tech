

--#########################################################
function ConsoleClass.ProcessEvents(this)
    this.consoleSizeNew=chiConsoleQueryUpdate();
    --print(this.consoleSizeNew,this.consoleSizeOld)
    if (this.consoleSizeNew>this.consoleSizeOld) then
        this.consoleSizeOld=this.consoleSizeNew;
        if (this.consoleSizeNew>0) then
            for k=1,this.consoleNumberOfDisplayLines do
                this.line[k].SetProperty(this.line[k],"Text",chiConsoleCopy(-k));
            end
             
        end
        
    end
end





--#########################################################
function ConsoleClass.SetProperty(this,property,value)
    if      (property=="Master") then
        this.master=value;
        
        value.slaveCount=value.slaveCount+1;
        k=value.slaveCount
        value.slaves[k]=this;
        
        this.SizeChanged(this);
    end
end

--#########################################################
function ConsoleClass.Redraw(this)
   
    for k=1,this.consoleNumberOfDisplayLines do
        --this.line[k].SetProperty(this.line[k],"Position",{this.xmin+this.paddingLeft,this.ymin+this.paddingBot+(k-1)*this.linespacing});
        this.line[k].SetProperty(this.line[k],"Position",{this.xmin+this.paddingLeft,this.ypos-this.paddingTop-(k)*this.linespacing});
        this.line[k].SetProperty(this.line[k],"Viewport",{this.xmin,this.ymin,this.xmax,this.ymax});
        this.line[k].SetProperty(this.line[k],"Text",chiConsoleCopy(-k));
    end
    
        
end





--#########################################################
function ConsoleClass.SizeChanged(this)
    
    if (this.master~=nil) then
        this.xmin=this.master.xmin+this.paddingLeft;
        this.ymin=this.master.ymin+this.paddingBot;
        this.xmax=this.master.xmax-this.paddingRight;
        this.ymax=this.master.ymax-this.paddingTop;
    end
    this.xpos=this.master.cursorX;
    this.ypos=this.master.cursorY;
    this.xSize=this.xmax-this.xmin+this.paddingLeft+this.paddingRight;
    --this.ySize=this.ymax-this.ymin+this.paddingBot+this.paddingTop-30;
    this.ySize=(this.consoleNumberOfDisplayLines+7)*this.linespacing
    --print(this.ySize)
    this.Redraw(this);
end


