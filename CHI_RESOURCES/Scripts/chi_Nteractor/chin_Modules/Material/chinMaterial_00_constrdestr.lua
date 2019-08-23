--######################################################### Constructor
function MaterialClass.New(toIndex)
    local this = setmetatable({},MaterialClass)
    MaterialCount = MaterialCount + 1; 
    
    this.index = toIndex;
    if (this.index==nil) then this.index=0; end
    chiMaterialQuery(5,this.index,this);
    
    --================================= Adding properties
    this.property    = {};
    this.property[0] = PropertyClass.New("Diffuse Color","Vec4");
    this.property[1] = PropertyClass.New("Ambient Color","Vec4");
    this.property[2] = PropertyClass.New("Specular Color","Vec4");
    this.property[3] = PropertyClass.New("Shininess","Scalar");
                    
    --================= Setting property values
    this.property[0].index = this.index;
    this.property[1].index = this.index;
    this.property[2].index = this.index;
    this.property[3].index = this.index;
    
    --================= Setting property values
    this.property[0].value = this.diffuse;
    this.property[1].value = this.ambient;
    this.property[2].value = this.specular;
    this.property[3].value = this.shininess;
    
    --================= Setting property display formats
    this.property[0].displayFormat = "%f,%f,%f,%f";
    this.property[1].displayFormat = "%f,%f,%f,%f";
    this.property[2].displayFormat = "%f,%f,%f,%f";
    this.property[3].displayFormat = "%.3f";
    
    --================= Setting property callback functions
    titem=this.property[0];
    function titem.CustomValueChanged(this)
        print("chiMaterialSetProperty",this.index,this.value[1],this.value[2],this.value[3],this.value[4])
        chiMaterialSetProperty(this.index,"Diffuse",this.value[1]/255,this.value[2]/255,this.value[3]/255,this.value[4]/255);
    end
    titem=this.property[1];
    function titem.CustomValueChanged(this)
        print("chiMaterialSetProperty")
        chiMaterialSetProperty(this.index,"Ambient",this.value[1]/255,this.value[2]/255,this.value[3]/255,this.value[4]/255);
    end
    titem=this.property[2];
    function titem.CustomValueChanged(this)
        print("chiMaterialSetProperty")
        chiMaterialSetProperty(this.index,"Specular",this.value[1]/255,this.value[2]/255,this.value[3]/255,this.value[4]/255);
    end
    titem=this.property[3];
    function titem.CustomValueChanged(this)
        print(this.value)
        print("chiMaterialSetProperty(",string.format(",%d,\"Shininess\",%d)",this.index,this.value))
        chiMaterialSetProperty(this.index,"Shininess",this.value);
    end
    
    return this;
end

