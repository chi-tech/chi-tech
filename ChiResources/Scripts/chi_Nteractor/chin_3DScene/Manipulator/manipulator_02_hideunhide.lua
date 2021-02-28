function ManipulatorClass.Hide(this)
    this.hidden = true;
    chiObjectSetProperty(this.xobjNum,"Hidden",true);
    chiObjectSetProperty(this.yobjNum,"Hidden",true);
    chiObjectSetProperty(this.zobjNum,"Hidden",true);
end

function ManipulatorClass.UnHide(this)
    this.hidden = false;
    chiObjectSetProperty(this.xobjNum,"Hidden",false);
    chiObjectSetProperty(this.yobjNum,"Hidden",false);
    chiObjectSetProperty(this.zobjNum,"Hidden",false);
end