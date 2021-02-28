/** \defgroup Lua_chil Y0 CHIL scripts 
 \ingroup LuaGeneralUtilities */
/** \defgroup Lua_Forms Forms
\ingroup Lua_chil
*/
/** Object for controlling forms.

A form object can exist as two different entities. It can either be
attached to a window, where it can be maximized/minimized as the user
wishes, or it can be internal to a window. The challenge of these
options is to cater for both needs.

\ingroup Lua_Forms
*/
class chilForm
{
public:
      char*               name;            ///< Form name (Will also be same as window)
      int                 sceneScope;      ///< Scene associated with form
      int                 displayerScope;  ///< Displayer associated with this form
      chilCameraObject    camera;          ///< Camera associated with this form
      chilForm            callbackReference; ///< Reference to self;
      callBackObject      callBackObject;  ///< Callback object for the form

      bool                initialized;     ///< Initialization flag
      int                 sceneScope;      ///< Scene this forms belongs to
      int                 displayerScope;  ///< Displayers this form belongs to

public:
      void        Initialize();
      void        callbackFunction();
};
/** Initializes basic form parameters.
\author Jan*/
void chilForm::Initialize(){}
/** Creates an external form.
\ingroup Lua_Forms
\author Jan*/
chilForm chilForms_CreateExternalForm(char name)
{
}
/** \brief Creates a basic orthographic camera for rendering to displayer coordinates.
This camera will follow the coordinates of the displayer its been created with.

\param name char Name of the new camera [optional]
\return chilCameraObject newCamera  Returns a camera object
\ingroup Lua_Camera
 -- */
void chilCreateOrthoDisplayerCamera(char name) {}
/** \defgroup Lua_Camera Camera control
\ingroup Lua_chil
*/

/** \brief Object to organize and register cameras.
Call example <I>chilCameraOrganizer.AddCamera(newCamera);</I>

\ingroup Lua_Camera
 -- */
class chilCameraOrganizer
{
public:
int itemCount;   ///< Number of cameras registered
void AddCamera(chilCameraObject newCamera);
};
/** \brief Call example <I>chilCameraOrganizer.AddCamera(newCamera);</I>

 Adds a camera to the organizer and registers the camera's callback.
 -- */
void chilCameraOrganizer::AddCamera(chilCameraObject newCamera){}
/** \brief Generic class for a camera. Override the <I>callbackFunction</I>
for customized response.

\ingroup Lua_Camera
 -- */
class chilCameraObject
{
public:
float         alpha      ;    ///< Camera orientation on the y axis
float         azimuth    ;    ///< Camera orientation on the x axis
float         tsi        ;    ///< Camera orientation on the x axis
float         x          ;
float         y          ;
float         z          ;
float         fov        ;

bool          left       ;
bool          right      ;
bool          up         ;
bool          down       ;
bool          forward    ;
bool          backward   ;
bool          lock1open  ;
bool          lock2open  ;
bool          zoomin     ;
bool          zoomout    ;

floar         zoomspeed  ;
floar         movespeed  ;
floar         orthoWidth ;

char          type;           ///< Either "Perspective" or "Orthographic"
float         xmax;
float         xmin;
float         ymax;
float         ymin;
float         zmax;
float         zmin;
int           sceneScope;
int           displayerScope;

public:
  void        Update(chilCameraObject self);
  void        callbackFunction(chilCameraObject self);
};
/** Updates camera type, position and rotation.
This function probably never needs to be changed.

\param chilCameraObject self Reference to the camera.

 -- */
void chilCameraObject::Update(chilCameraObject self) {}
/** Function that controls the behavior of the camera.

\param chilCameraObject self Reference to the camera.

 -- */
void chilCameraObject::callbackFunction(chilCameraObject self) {}
/** \brief Creates a basic orthographic camera for rendering to window coordinates.
This camera will follow the coordinates of the scene/window its been created in.

\param name char Name of the new camera [optional]
\return chilCameraObject newCamera  Returns a camera object
\ingroup Lua_Camera
 -- */
void chilCreateOrthoWindowCamera(char name) {}
/** \brief Creates a basic thirdperson camera.

\param name char Name of the new camera [optional]
\return chilCameraObject newCamera  Returns a camera object
\ingroup Lua_Camera
 -- */
void chilCreateThirdPersonCamera(char name) {}
/** \brief Creates a camera that revolves around a point.

\param name char Name of the new camera [optional]
\return chilCameraObject newCamera  Returns a camera object
\ingroup Lua_Camera
 -- */
void chilCreateRevolverCamera(char name) {}
/** \defgroup Lua_Callbacks Callbacks
\ingroup Lua_chil

Scripted object for handling **callbacks**. Requires that an object has function
called **callbackFunction**.

 ## Example
\code
SomeObject = {};
SomeObject.callbackFunction = function(this)
  --Handle events
end
--
callBackObj = chilCallbacks.MakeCallback(SomeObject);
chilCallbacks.PushCallback(callBackObj);
\endcode
*/

/** Structure for defining a callback object

 */
struct callBackObject
{
      chilObject  parent;
      bool        enabled;        ///< Enabled/Disabled flag
      bool        once;           ///< Flag for terminating callback after first call [default: false]
      float       cyclic;         ///< Cyclical execution [default: nil]
};

/** Contains all the functionality to implement callbacks.

\ingroup Lua_Callbacks
 */
class chilCallbacks
{
public:
      int count; ///<Default 0

      void                PushCallback(callBackObject callBack);
      callBackObject      MakeCallback(chilObject parentObject);
      void                Execute();
};

/** Adds a callback object to the callback stack.

\param callBack The object to be pushed.
 */
void chiCallbacks::PushCallback(callBackObject callBack){}

/** Creates a callBack object for the given parent object.
The parent object needs to have a method called "callbackFunction"
and two properties: callbackReference and callbackParameters

\param     parentObject The object to be pushed.
\return    callBackObject Returns a callBackObject.
 */
callBackObject chiCallbacks::MakeCallback(chilObject parentObject){}

/** Executes all callback functions

 */
void chiCallbacks::Execute(){}

};
