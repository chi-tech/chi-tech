--========================================================= Selection class
SelectionClass = {}
SelectionClass.__index = SelectionClass;

SELECTION_TYPE_VOID             = 0;
SELECTION_TYPE_OBJECT           = 1;
SELECTION_TYPE_MATERIAL         = 2;
SELECTION_TYPE_TRANSFORMATION   = 3;
SELECTION_TYPE_PROPERTY         = 4;
SELECTION_TYPE_TEXTURE          = 5;
SELECTION_TYPE_SHADER           = 6;
SELECTION_TYPE_TABLE            = 7;

SELECTION_ORIGIN_3D       = 0;
SELECTION_ORIGIN_TREEVIEW = 1;
SELECTION_ORIGIN_CODE     = 2;

function SelectionClass.New()
    local this=setmetatable({},SelectionClass)
        
    this.type     =SELECTION_TYPE_VOID;
    this.index    = -1;
    this.item     = nil;
    this.origin   = SELECTION_ORIGIN_CODE;
    this.originFeature = nil;
    
    return this;
end