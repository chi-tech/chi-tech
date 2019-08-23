--#{
--/** \defgroup Lua_Camera Camera control
--\ingroup Lua_chil
--*/
--
--#}


--#{
--/** \brief Object to organize and register cameras.
--Call example <I>chilCameraOrganizer.AddCamera(newCamera);</I>
--
--\ingroup Lua_Camera
-- -- */
--class chilCameraOrganizer
--{
--public:
--int itemCount;   ///< Number of cameras registered
--void AddCamera(chilCameraObject newCamera);
--};
--#}
chilCameraOrganizer = {}
chilCameraOrganizer.itemCount = 0;


--#{
--/** \brief Call example <I>chilCameraOrganizer.AddCamera(newCamera);</I>
--
-- Adds a camera to the organizer and registers the camera's callback.
-- -- */
--void chilCameraOrganizer::AddCamera(chilCameraObject newCamera){}
--#}
chilCameraOrganizer.AddCamera = function (cameraObj)
    chilCameraOrganizer.itemCount = chilCameraOrganizer.itemCount+1;
    local index = chilCameraOrganizer.itemCount;
    chilCameraOrganizer[index] = cameraObj;

    local newCallback = chilCallbacks.MakeCallback(cameraObj);
    chilCallbacks.PushCallback(newCallback);
end

--#{
--/** \brief Generic class for a camera. Override the <I>callbackFunction</I>
--for customized response.
--
--\ingroup Lua_Camera
-- -- */
--class chilCameraObject
--{
--public:
--float         alpha      ;    ///< Camera orientation on the y axis
--float         azimuth    ;    ///< Camera orientation on the x axis
--float         tsi        ;    ///< Camera orientation on the x axis
--float         x          ;
--float         y          ;
--float         z          ;
--float         fov        ;
--
--bool          left       ;
--bool          right      ;
--bool          up         ;
--bool          down       ;
--bool          forward    ;
--bool          backward   ;
--bool          lock1open  ;
--bool          lock2open  ;
--bool          zoomin     ;
--bool          zoomout    ;
--
--floar         zoomspeed  ;
--floar         movespeed  ;
--floar         orthoWidth ;
--
--char          type;           ///< Either "Perspective" or "Orthographic"
--float         xmax;
--float         xmin;
--float         ymax;
--float         ymin;
--float         zmax;
--float         zmin;
--int           sceneScope;
--int           displayerScope;
--
--public:
--  void        Update(chilCameraObject self);
--  void        callbackFunction(chilCameraObject self);
--};
--#}

--#{
--/** Updates camera type, position and rotation.
--This function probably never needs to be changed.
--
--\param chilCameraObject self Reference to the camera.
--
-- -- */
--void chilCameraObject::Update(chilCameraObject self) {}
--#}

--#{
--/** Function that controls the behavior of the camera.
--
--\param chilCameraObject self Reference to the camera.
--
-- -- */
--void chilCameraObject::callbackFunction(chilCameraObject self) {}
--#}