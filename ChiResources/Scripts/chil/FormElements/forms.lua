




--#{
--/** \defgroup Lua_Forms Forms
--\ingroup Lua_chil
--*/
--#}

--#{
--/** Object for controlling forms.
--
--A form object can exist as two different entities. It can either be
--attached to a window, where it can be maximized/minimized as the user
--wishes, or it can be internal to a window. The challenge of these
--options is to cater for both needs.
--
--\ingroup Lua_Forms
--*/
--class chilForm
--{
--public:
--      char*               name;            ///< Form name (Will also be same as window)
--      int                 sceneScope;      ///< Scene associated with form
--      int                 displayerScope;  ///< Displayer associated with this form
--      chilCameraObject    camera;          ///< Camera associated with this form
--      chilForm            callbackReference; ///< Reference to self;
--      callBackObject      callBackObject;  ///< Callback object for the form
--
--      bool                initialized;     ///< Initialization flag
--      int                 sceneScope;      ///< Scene this forms belongs to
--      int                 displayerScope;  ///< Displayers this form belongs to
--
--public:
--      void        Initialize();
--      void        callbackFunction();
--};
--#}

--#{
--/** Initializes basic form parameters.
--\author Jan*/
--void chilForm::Initialize(){}
--#}


--#{
--/** Creates an external form.
--\ingroup Lua_Forms
--\author Jan*/
--chilForm chilForms_CreateExternalForm(char name)
--{
--}
--#}
function chilForms_CreateExternalForm(name)
    local newForm = {}
    newForm.name                = name;
    newForm.sceneScope          = -1;
    newForm.displayerScope      = -1;
    newForm.camera              = -1;
    newForm.callbackReference   = newForm;
    newForm.callBackObject      = nil;

    newForm.initialized         = false;

    --====================================================== Set initial values for displayer
    newForm.sceneScope          = chiWindowCreate(name, true);
    newForm.displayerScope      = 0;

    --====================================================== Initialize
    newForm.Initialize = function (this)
        curScene, curDispl = chiGetScene();
        chiBindScene(this.sceneScope,this.displayerScope);
        this.camera  = chilCreateOrthoWindowCamera(this.name.." camera");
        chilCameraOrganizer.AddCamera(this.camera);



        local text1 = chiTextCreate("Text");
        chiTextSetProperty(text1,1,"Hello World");

        chiBindScene(curScene,curDispl);
        this.initialized=true;
    end

    --====================================================== Callback function
    newForm.callbackFunction = function (this)
        if ((this.initialized~=true) and (chiGetSceneCount()>=(this.sceneScope + 1))) then
            this.Initialize(this);
        end

    end

    --====================================================== Callback object
    newForm.callBackObject = chilCallbacks.MakeCallback(newForm);
    chilCallbacks.PushCallback(newForm.callBackObject);

    return newForm;
end