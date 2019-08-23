



function chilForms_CreateVSplit(name)
    local newVSplit = {}
    newVSplit.name                = name;
    newVSplit.sceneScope          = -1;
    newVSplit.displayerScope      = -1;
    --.camera
    newVSplit.callbackReference   = newVSplit;
    newVSplit.callBackObject      = nil;

    --newVSplit.initialized         = false;

    newVSplit.anchors = {left=nil,right=nil,top=nil,bottom=nil }
    newVSplit.parent  = nil;
    newVSplit.children = {itemCount=0 }
    newVSplit.selected = false;

    newVSplit.position = {x=100,y=100};

    newVSplit.lastObject          = -1;

    newVSplit.sceneScope, newVSplit.displayerScope = chiGetScene();

    newVSplit.divider = chi3DLineCreate("");

    chi3DLineAddVertex(newVSplit.divider,100.0,100.0, 10.0);
    chi3DLineAddVertex(newVSplit.divider,100.0,400.0, 10.0);
    chi3DLineSetStipple(newVSplit.divider,false,1,5.0);

    --====================================================== Callback function
    newVSplit.callbackFunction = function (this)
        --=================================== Mouse move event
        if (WM_MOUSEMOVE.occured) then
            curScene, curDispl = chiGetScene();
            chiBindScene(this.sceneScope,this.displayerScope);

            --========================= Pointer for divider
            if ((WM_MOUSEMOVE.iPar5 == this.divider)) then
                chiSetWindowProperties("CURSOR",108);
            end
            if ((WM_MOUSEMOVE.iPar5 ~= this.divider) and (lastObject == this.divider)) then
                chiSetWindowProperties("CURSOR",-1);
            end
            if (this.selected == true) then
                this.position.x = this.position.x + (WM_MOUSEMOVE.iPar0 - WM_MOUSEMOVE.iPar2);
                this.position.y = this.position.y + (WM_MOUSEMOVE.iPar1 - WM_MOUSEMOVE.iPar3);

                if (this.anchors.top   ==nil) then chi3DLineChangeVertex(this.divider,1,this.position.x,400.0,10.0); end
                if (this.anchors.bottom==nil) then chi3DLineChangeVertex(this.divider,0,this.position.x,  0.0,10.0); end
            end

            lastObject = WM_MOUSEMOVE.iPar5;
            chiBindScene(curScene,curDispl);
        end

        --=================================== Mouse down event
        if (WM_LBUTTONDOWN.occured and (WM_LBUTTONDOWN.iPar5 == this.divider)) then
            chilSelectionStack.PushSelection(this);
        end

        --=================================== Mouse up event
        if (WM_LBUTTONUP.occured) then
            chilSelectionStack.DeselectItem(this);
        end


        --=================================== Mouse move event
        if (WM_SIZE.occured) then
            if (this.anchors.top   ==nil) then chi3DLineChangeVertex(this.divider,1,this.position.x,400.0,10.0); end
            if (this.anchors.bottom==nil) then chi3DLineChangeVertex(this.divider,0,this.position.x,  0.0,10.0); end
        end
    end

    --====================================================== Callback object
    newVSplit.callBackObject = chilCallbacks.MakeCallback(newVSplit);
    chilCallbacks.PushCallback(newVSplit.callBackObject);

    --====================================================== Select method
    newVSplit.Select = function (this)
        this.selected = true;
    end

    --====================================================== DeSelect method
    newVSplit.DeSelect = function (this)
        this.selected = false;
    end

    return newVSplit;
end