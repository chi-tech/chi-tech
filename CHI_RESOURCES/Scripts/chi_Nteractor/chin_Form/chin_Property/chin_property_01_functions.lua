--######################################################### Value Changed
function PropertyClass.ValueChanged(this)
    if (this.CustomValueChanged~=nil) then
        this.CustomValueChanged(this);
    end
end
