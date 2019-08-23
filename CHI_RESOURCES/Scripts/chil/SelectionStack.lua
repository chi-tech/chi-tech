





chilSelectionStack = {itemCount=0}

chilSelectionStack.PushSelection = function (object)
    chilSelectionStack.itemCount = chilSelectionStack.itemCount + 1;
    local index = chilSelectionStack.itemCount;
    chilSelectionStack[index] = object;
    object.Select(object);
end

chilSelectionStack.Clear = function ()
    for k=1,chilSelectionStack.itemCount do
        if (chilSelectionStack[k] ~= nil) then
            if (chilSelectionStack[k].DeSelect ~= nil) then
                chilSelectionStack[k].DeSelect(chilSelectionStack[k]);
            end
            chilSelectionStack[k]=nil;
        end
    end
    chilSelectionStack.itemCount = 0;
end

chilSelectionStack.DeselectItem = function (object)
    for k=1,chilSelectionStack.itemCount do
        if (chilSelectionStack[k] ~= nil) then
            if (chilSelectionStack[k] == object) then
                if (chilSelectionStack[k].DeSelect ~= nil) then
                    chilSelectionStack[k].DeSelect(chilSelectionStack[k]);
                end
                chilSelectionStack[k]=nil;
            end
        end
    end
end