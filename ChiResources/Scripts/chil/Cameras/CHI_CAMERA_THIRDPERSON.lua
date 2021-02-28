--######################################################### Camera def
--#{
--/** \brief Creates a basic thirdperson camera.
--
--\param name char Name of the new camera [optional]
--\return chilCameraObject newCamera  Returns a camera object
--\ingroup Lua_Camera
-- -- */
--void chilCreateThirdPersonCamera(char name) {}
--#}
function chilCreateThirdPersonCamera(name)
    local newCamera = {}
    newCamera.alpha        = 0;        -- camera orientation on the x axis
    newCamera.azimuth      = 0;    -- camera orientation on the y axis
    newCamera.tsi          = 0;        -- something
    newCamera.x            = 0.0;      -- x pos. 
    newCamera.y            = 0.0;    -- y pos.
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
    newCamera.orthoWidth   = 50.0;
    newCamera.type         = "Perspective"
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
        else
            chiGraphicsCameraType(1);
        end
        chiGraphicsOrientCamera(0,this.alpha,this.azimuth,this.tsi);
        chiGraphicsPositionCamera(this.x,this.y,this.z);

        chiBindScene(currentScene,currentDisplayer);
    end
    
    --===================================================== Callback function
    newCamera.callbackFunction = function(this)
        --print(chi_frameRate);
        if (WM_KEYDN.occured and (WM_KEYDN.iPar4 == this.sceneScope)) then
            --print(string.format("WM_KEYDN %6.4f",chi_programTime));
            if (WM_KEYDN.iPar0 == 65) then this.left = true; end
            if (WM_KEYDN.iPar0 == 68) then this.right = true; end
            if (WM_KEYDN.iPar0 == 87) then this.forward = true; end
            if (WM_KEYDN.iPar0 == 83) then this.backward = true; end
            if (WM_KEYDN.iPar0 == 69) then this.up = true; end
            if (WM_KEYDN.iPar0 == 81) then this.down = true; end
        end

        if (WM_MOUSEMOVE.occured) then
            --print(string.format("Object selected: %d",WM_MOUSEMOVE.iPar5));
        end

    
        if (WM_KEYUP.occured or WM_MOUSELEAVE.occured) then
            --print(string.format("WM_KEYUP %6.4f",chi_programTime));
            if (WM_KEYUP.iPar0 == 65) then this.left = false; end
            if (WM_KEYUP.iPar0 == 68) then this.right = false; end
            if (WM_KEYUP.iPar0 == 87) then this.forward = false; end
            if (WM_KEYUP.iPar0 == 83) then this.backward = false; end
            if (WM_KEYUP.iPar0 == 69) then this.up = false; end
            if (WM_KEYUP.iPar0 == 81) then this.down = false; end
            
            if (WM_MOUSELEAVE.occured) then 
                    this.left = false; 
                    this.right = false;
                    this.forward = false;
                    this.backward = false;
                    this.up = false;
                    this.down = false;
                    this.lock2open=false;
            end
        end
      
        
        local dx = 0; 
        local dy = 0;
        local dz = 0;
        
        
        if ((this.left)    ) then dx = 0.0 - this.movespeed; end
        if ((this.right)   ) then dx = 0.0 + this.movespeed;  end
        if ((this.forward) ) then dy = 0.0 + this.movespeed;  end
        if ((this.backward)) then dy = 0.0 - this.movespeed;  end
        --if ((this.up)      ) then   end
        --if ((this.down)    ) then   end
        --
        this.x = this.x + dx;
        this.y = this.y + dy;
        this.z = this.z + dz;
        
        if (this.x>this.xmax) then this.x = this.x + dx; end
        if (this.x<this.xmin) then this.x = this.x + dx; end
        if (this.y>this.ymax) then this.y = this.y + dy; end
        if (this.y<this.ymin) then this.y = this.y + dy; end
        if (this.z>this.zmax) then this.z = this.z + dz; end
        if (this.z<this.zmin) then this.z = this.z + dz; end
        
        this.Update(this);
    end
    newCamera.CameraControl = newCamera.callbackFunction
    
    print("Third person camera made")
    return newCamera;
end


