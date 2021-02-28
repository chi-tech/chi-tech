--######################################################### Constructor
function PropertyClass.New(name,propType)
    local this=setmetatable({},PropertyClass);
    
    
    
    if     (propType=="Name") then
        this.name=name; this.propType=1;
    elseif (propType=="Color") then
        this.name=name; this.propType=2; 
        this.r=0.0;
        this.g=0.0;
        this.b=0.0;
        this.a=1.0;
    elseif (propType=="Scalar") then
        this.name=name; this.propType=3;
        this.value=0.0;
        this.displayFormat="%f";
    elseif (propType=="Vec2") then
        this.name=name; this.propType=4;
        this.value={0.0,0.0};
        this.displayFormat="%f,%f";
    elseif (propType=="Vec3") then
        this.name=name; this.propType=5;
        this.value={0.0,0.0,0.0};
        this.displayFormat="%f,%f,%f";
    elseif (propType=="Vec4") then
        this.name=name; this.propType=6;
        this.value={0.0,0.0,0.0,0.0};
        this.displayFormat="%f,%f,%f,%f";
    elseif (propType=="Droplist") then
        this.name=name; this.propType=7;
        this.value={}
        this.valueCount = 0;
        this.displayFormat = "%f"
        this.initOption = 0;
        this.option = 0;
    elseif (propType=="Boolean") then
        this.name=name; this.propType=8;
        this.value = false;
    end
    
    return this;
end