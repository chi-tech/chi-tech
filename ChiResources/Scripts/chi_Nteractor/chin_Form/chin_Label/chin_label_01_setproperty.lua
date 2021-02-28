--#########################################################
function LabelClass.SetProperty(this,property,value)
    if     (property=="Text") then
        this.text=value;
        r=this.fontColor[1];
        g=this.fontColor[2];
        b=this.fontColor[3];
        a=this.fontColor[4];
        chiSetLabel3D(this.textNum,this.text,this.xpos+this.paddingLeft,this.ypos+this.paddingBot,r,g,b,a,this.fontType)
    elseif (property=="Position") then
        if (this.float) then
            this.xpos=value[1];
            this.ypos=value[2];
            chiSetLabel3D(this.textNum,this.text,this.xpos+this.paddingLeft,this.ypos+this.paddingBot,r,g,b,a,this.fontType)
        end
    elseif (property=="Viewport") then
        xmin=value[1];
        ymin=value[2];
        xmax=value[3];
        ymax=value[4];
        chiSetLabelProperty3D(this.textNum,"Viewport",xmin,ymin,xmax,ymax);
        chiViewportSetProperty(this.obj1View,xmin,ymin,xmax,ymax);
    elseif (property=="ViewportEnable") then
        chiSetLabelProperty3D(this.textNum,"ViewportEnable",value);
    elseif (property=="Master") then
        this.master=value;
        
        value.slaveCount=value.slaveCount+1;
        k=value.slaveCount
        value.slaves[k]=this;
     elseif (property=="Color") then
        this.fontColor[1]=value[1];
        this.fontColor[2]=value[2];
        this.fontColor[3]=value[3];
        this.fontColor[4]=value[4];
        r=value[1];
        g=value[2];
        b=value[3];
        a=value[4];
        chiSetLabel3D(this.textNum,this.text,this.xpos+this.paddingLeft,this.ypos+this.paddingBot,r,g,b,a,this.fontType)
     elseif (property=="Float") then
        this.float=value;
    end
end