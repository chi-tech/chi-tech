function TableClass.New(name)
    local this = setmetatable({},TableClass);
    
    local objName="Table"..string.format("%03d",TableCount);
    this.name = name;
    
    this.tableIndex = chiTableCreate(objName);
    
    return this;
end