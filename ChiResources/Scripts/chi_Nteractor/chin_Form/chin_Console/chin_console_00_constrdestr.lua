

--#########################################################
function ConsoleClass.New(name)
    local this=setmetatable({},ConsoleClass)
    ConsoleCount=ConsoleCount+1
    
    objName="Console"..string.format("%02d",ConsoleCount);
    
    this.name=name;
    this.xpos=100;
    this.ypos=100;
    this.xmin=100;
    this.ymin=100;
    this.xmax=200;
    this.ymax=200;
    this.ySize=200;
    this.xSize=200;
    this.linespacing=15;
    this.paddingTop=10;
    this.paddingBot=0;
    this.paddingLeft=5;
    this.paddingRight=5;
    this.consoleSizeNew=0;
    this.consoleSizeOld=0;
    this.consoleNumberOfDisplayLines=100;
    --this.consoleInitialized=false;
    
    this.line={}
    for k=1,this.consoleNumberOfDisplayLines do
        this.line[k]={}
        this.line[k]=LabelClass.New(objName.."Label001");
        this.line[k].fontType=2;
        --this.line[k].SetProperty(this.line[k],"Position",{this.xmin+this.paddingLeft,this.ymin+this.paddingBot+(k-1)*this.linespacing});
        this.line[k].SetProperty(this.line[k],"Position",{this.xmin+this.paddingLeft,this.ymax-this.paddingTop-(k-1)*this.linespacing});
        this.line[k].SetProperty(this.line[k],"Viewport",{this.xmin,this.ymin,this.xmax,this.ymax});
        this.line[k].SetProperty(this.line[k],"ViewportEnable",true);
        this.line[k].SetProperty(this.line[k],"Text",chiConsoleCopy(-k));
    end
    
    
    return this;
end