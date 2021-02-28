--######################################################### Constructor
function TextureClass.New(toIndex)
    local this = setmetatable({},TextureClass)
    TextureCount = TextureCount + 1; 
    
    this.index = toIndex;
    if (this.index==nil) then this.index=0; end
    chiTextureQuery(5,this.index,this);
    
    --================================= Adding properties
    this.property    = {};
    this.property[0] = PropertyClass.New("Diffuse Color","Name");
    
    --================= Setting property values
    this.property[0].index = this.index;
    
    --================= Setting property values
    this.property[0].name = this.name;
    
    --================= Setting property callback functions
    titem=this.property[0];
    function titem.CustomValueChanged(this)
       print("No script set");
    end
    
    
    return this;
end