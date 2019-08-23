--######################################################### Camera def
--#{
--/** \brief Creates a basic orthographic camera for rendering to window coordinates.
--This camera will follow the coordinates of the scene/window its been created in.
--
--\param name char Name of the new camera [optional]
--\return chilCameraObject newCamera  Returns a camera object
--\ingroup Lua_Camera
-- -- */
--void chilCreateOrthoWindowCamera(char name) {}
--#}
function chilCreateOrthoWindowCamera(name)
    local newCamera = {}
    newCamera.alpha        = 0;        -- camera orientation on the x axis
    newCamera.azimuth      = 0;        -- camera orientation on the y axis
    newCamera.tsi          = 0;        -- something
    newCamera.x            = 0.0;      -- x pos.
    newCamera.y            = 0.0;     -- y pos.
    newCamera.z            = 10.0;      -- z pos.
    newCamera.fov          = 40.0;
    newCamera.left         = false;
    newCamera.right        = false;
    newCamera.up           = false;
    newCamera.down         = false;
    newCamera.forward      = false;
    newCamera.backward     = false;
    newCamera.lock1open    = false;
    newCamera.lock2open    = false;
    newCamera.zoomin       = false;
    newCamera.zoomout      = false;
    newCamera.zoomspeed    = 2.0;
    newCamera.movespeed    = 1.0*0.0166666;
    newCamera.orthoWidth   = 200.0;
    newCamera.type         = "Orthographic"
    newCamera.xmax         = 10000;
    newCamera.xmin         =-10000;
    newCamera.ymax         = 10000;
    newCamera.ymin         =-10000;
    newCamera.zmax         = 10000;
    newCamera.zmin         =-10000;
    newCamera.sceneScope,newCamera.displayerScope   = chiGetScene();

    newCamera.callbackReference = newCamera;

    --===================================================== Update function
    newCamera.Update = function(this)
        currentScene,currentDisplayer = chiGetScene();
        chiBindScene(this.sceneScope,this.displayerScope);

        if (this.type=="Orthographic") then
            chiGraphicsCameraType(2);
            chiGraphicsCameraOrthoWidth(this.orthoWidth,true);
        else
            chiGraphicsCameraType(1);
        end
        chiGraphicsOrientCamera(0,this.alpha,this.azimuth,this.tsi);
        chiGraphicsPositionCamera(this.x,this.y,this.z);

        chiBindScene(currentScene,currentDisplayer);
    end

    --===================================================== Callback function
    newCamera.callbackFunction = function(this)

        if (WM_SIZE.occured) then
            currentScene,currentDisplayer = chiGetScene();
            chiBindScene(this.sceneScope,this.displayerScope);
            xSize,ySize,x,y =chiGetWindowProperties();
            this.orthoWidth = xSize;
            chiBindScene(currentScene,currentDisplayer);
            this.Update(this);
        end


    end
    newCamera.CameraControl = newCamera.callbackFunction

    print("Ortho window camera made")
    return newCamera;
end

