--#{
--/** \defgroup Lua_Callbacks Callbacks
--\ingroup Lua_chil
--
--Scripted object for handling **callbacks**. Requires that an object has function
--called **callbackFunction**.
--
-- ## Example
--\code
--SomeObject = {};
--SomeObject.callbackFunction = function(this)
--  --Handle events
--end
----
--callBackObj = chilCallbacks.MakeCallback(SomeObject);
--chilCallbacks.PushCallback(callBackObj);
--\endcode
--*/
--
--/** Structure for defining a callback object
--
-- */
--struct callBackObject
--{
--      chilObject  parent;
--      bool        enabled;        ///< Enabled/Disabled flag
--      bool        once;           ///< Flag for terminating callback after first call [default: false]
--      float       cyclic;         ///< Cyclical execution [default: nil]
--};
--
--#}

--#{
--/** Contains all the functionality to implement callbacks.
--
--\ingroup Lua_Callbacks
-- */
--class chilCallbacks
--{
--public:
--      int count; ///<Default 0
--
--      void                PushCallback(callBackObject callBack);
--      callBackObject      MakeCallback(chilObject parentObject);
--      void                Execute();
--};
--
--#}

chilCallbacks={}
chilCallbacks.count = 0;

--=============================================== Pushes a Callback to the stack
--#{
--/** Adds a callback object to the callback stack.
--
--\param callBack The object to be pushed.
-- */
--void chiCallbacks::PushCallback(callBackObject callBack){}
--
--#}
chilCallbacks.PushCallback = function (value)
    chilCallbacks.count = chilCallbacks.count+1;
    local index = chilCallbacks.count
    chilCallbacks[index] = value;
end

--=============================================== Creates a standard Callback
--#{
--/** Creates a callBack object for the given parent object.
--The parent object needs to have a method called "callbackFunction"
--and two properties: callbackReference and callbackParameters
--
--\param     parentObject The object to be pushed.
--\return    callBackObject Returns a callBackObject.
-- */
--callBackObject chiCallbacks::MakeCallback(chilObject parentObject){}
--
--#}
chilCallbacks.MakeCallback = function (parentObject)
    local newCallBack = {}
    newCallBack.parent   = parentObject;
    newCallBack.enabled  = true;
    newCallBack.once     = false;
    newCallBack.cyclic   = nil;
    
    return newCallBack;
end

--=============================================== Execute callbacks
--#{
--/** Executes all callback functions
--
-- */
--void chiCallbacks::Execute(){}
--
--#}
chilCallbacks.Execute = function ()
    for k=1,chilCallbacks.count,1 do
        local callbackObject = chilCallbacks[k].parent;
        if (not (callbackObject==nil)) then
            
            if (not (callbackObject.callbackFunction==nil)) then
                if (callbackObject.callbackReference==nil) then
                    callbackObject.callbackFunction();
                else
                    if (callbackObject.callbackParameters==nil) then
                        callbackObject.callbackFunction(callbackObject.callbackReference);
                    else
                        local ref    = callbackObject.callbackReference;
                        local params = callbackObject.callbackParameters;
                        callbackObject.callbackFunction(ref, params);
                    end
                end
                
            end
        end
    end
end

