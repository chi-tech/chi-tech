--========================================================= Selection Stack Management
selectionStack={}
selectionStack.itemCount=0;
selectionStack.__index = selectionStack;
selectionStack.PushCallBacks={};
selectionStack.PushCallBackCount=0;
selectionStack.PullCallBacks={};
selectionStack.PullCallBackCount=0;
selectionStack.CtlDown   = false;
selectionStack.ShiftDown = false;

--========================================================= Push item to the stack
function selectionStack.PushItem(selectionItem)
    if (selectionStack.CtlDown) then
        selectionStack.itemCount=selectionStack.itemCount+1;
        selectionStack[selectionStack.itemCount]=selectionItem;

        --=================== Call registered callbacks
        for k=1,selectionStack.PushCallBackCount do
            selectionStack.PushCallBacks[k](selectionItem);
        end
    else
        selectionStack.Clear()
        selectionStack.itemCount=selectionStack.itemCount+1;
        selectionStack[selectionStack.itemCount]=selectionItem;

        --=================== Call registered callbacks
        for k=1,selectionStack.PushCallBackCount do
            selectionStack.PushCallBacks[k](selectionItem);
        end
    end
    return selectionStack.itemCount;
end

--========================================================= Pull item form the stack
function selectionStack.PullItem(selIndex)
    
end

--========================================================= Clear the stack
function selectionStack.Clear()
    for m=1,selectionStack.itemCount do
        selectionItem=selectionStack[m];
        for k=1,selectionStack.PullCallBackCount do
            selectionStack.PullCallBacks[k](selectionItem);
        end
    end
    selectionStack.itemCount=0;
end

--========================================================= Push Selection Callback
function selectionStack.PushSelectionCallback(functionReference)
    selectionStack.PushCallBackCount=selectionStack.PushCallBackCount+1;
    k=selectionStack.PushCallBackCount;
    selectionStack.PushCallBacks[k]=functionReference;
end

function selectionStack.PushDeSelectionCallback(functionReference)
    selectionStack.PullCallBackCount=selectionStack.PullCallBackCount+1;
    k=selectionStack.PullCallBackCount;
    selectionStack.PullCallBacks[k]=functionReference;
end

--========================================================= Selection stack keyboard events
function selectionStack.ProcessEvents()
    if (WM_CTLBUTTONDN.occured) then selectionStack.CtlDown = true; end
    if (WM_CTLBUTTONUP.occured) then selectionStack.CtlDown = false; end
end

dRegisterFeature(selectionStack);

dofile(chinSelectionStackDir.."chin_selectionclass.lua");
