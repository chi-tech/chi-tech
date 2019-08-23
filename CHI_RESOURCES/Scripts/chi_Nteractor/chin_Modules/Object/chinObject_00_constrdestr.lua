--######################################################### Constructor
function ObjectClass.New(toIndex)
    local this = setmetatable({},ObjectClass)
    ObjectCount = ObjectCount + 1; 
    
    this.index = toIndex;
    if (this.index==nil) then this.index=0; end
    chiObjectQuery(5,this.index,this);
    
    

    
    return this;
end