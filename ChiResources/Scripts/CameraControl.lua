Camera_alpha        = 0;        -- camera orientation on the x axis
Camera_azimuth      = -90.0;    -- camera orientation on the y axis
Camera_tsi          = 0;        -- something
Camera_x            = 0.0;      -- x pos. 
Camera_y            = -10.0;    -- y pos.
Camera_z            = 0.0;      -- z pos.
Camera_fov          = 40.0;
Camera_left         = false;
Camera_right        = false;
Camera_up           = false;
Camera_down         = false;
Camera_forward      = false; 
Camera_backward     = false;
Camera_lock1open    = false;
Camera_lock2open    = false;
Camera_zoomin       = false;
Camera_zoomout      = false;
Camera_zoomspeed    = 2.0;
Camera_movespeed    = 1.0;
Camera_orthoWidth   = 50.0;
Camera_type         = "Perspective";
Camera_xmax         = 10000;
Camera_xmin         =-10000;
Camera_ymax         = 10000;
Camera_ymin         =-10000;
Camera_zmax         = 10000;
Camera_zmin         =-10000;


chiGraphicsOrientCamera(0,Camera_alpha,Camera_azimuth,Camera_tsi);
chiGraphicsPositionCamera(Camera_x,Camera_y,Camera_z);

function CameraControl()

    if (WM_KEYDN.occured) then
       if (WM_KEYDN.iPar0 == 65) then Camera_left = true; end
       if (WM_KEYDN.iPar0 == 68) then Camera_right = true; end
       if (WM_KEYDN.iPar0 == 87) then Camera_forward = true; end
       if (WM_KEYDN.iPar0 == 83) then Camera_backward = true; end
       if (WM_KEYDN.iPar0 == 69) then Camera_up = true; end
       if (WM_KEYDN.iPar0 == 81) then Camera_down = true; end
    end

    if (WM_KEYUP.occured or WM_MOUSELEAVE.occured) then
        --print("WM_KEYUP");
       if (WM_KEYUP.iPar0 == 65) then Camera_left = false; end
       if (WM_KEYUP.iPar0 == 68) then Camera_right = false; end
       if (WM_KEYUP.iPar0 == 87) then Camera_forward = false; end
       if (WM_KEYUP.iPar0 == 83) then Camera_backward = false; end
       if (WM_KEYUP.iPar0 == 69) then Camera_up = false; end
       if (WM_KEYUP.iPar0 == 81) then Camera_down = false; end
       
       if (WM_MOUSELEAVE.occured) then 
            Camera_left = false; 
            Camera_right = false;
            Camera_forward = false;
            Camera_backward = false;
            Camera_up = false;
            Camera_down = false;
            Camera_lock2open=false;
            print("Muff");
       end
    end
    
    if (WM_KEYDN.occured and (WM_KEYDN.iPar0 == 84)) then
        if (Camera_type=="Perspective") then
            Camera_type="Orthographic";
            chiGraphicsCameraType(2);
        else
            Camera_type="Perspective";
            chiGraphicsCameraType(1);
        end
    end
    
    if (WM_CTLBUTTONDN.occured) then Camera_lock1open = true; end
    --if ((WM_LBUTTONDOWN.occured) or (WM_MBUTTONDOWN.occured)) then Camera_lock2open = true; end
    if ((WM_LBUTTONDOWN.occured) ) then Camera_lock2open = true; end
    
    if (WM_CTLBUTTONUP.occured) then Camera_lock1open = false; end
    --if ((WM_LBUTTONUP.occured) or (WM_MBUTTONUP.occured)) then Camera_lock2open = false; end
    if ((WM_LBUTTONUP.occured) ) then Camera_lock2open = false; end
    
    if ((true) and (Camera_lock2open) and (WM_MOUSEMOVE.occured)) then
        mouse_dx = WM_MOUSEMOVE.iPar0 - WM_MOUSEMOVE.iPar2;
        mouse_dy = WM_MOUSEMOVE.iPar1 - WM_MOUSEMOVE.iPar3;
        Camera_alpha=Camera_alpha-0.2*mouse_dx;
        Camera_azimuth=Camera_azimuth-0.2*mouse_dy;
        chiGraphicsOrientCamera(0,Camera_alpha,Camera_azimuth,Camera_tsi);
    end
    
    --if (WM_MBUTTONDOWN.occured) then
    --    Camera_zoomin       = true;
    --    Camera_zoomout      = false;
    --end
    --
    --if (WM_MBUTTONUP.occured) then
    --    Camera_zoomin       = false;
    --    Camera_zoomout      = true;
    --end
    
    if (Camera_zoomin) then
        Camera_fov=Camera_fov-Camera_zoomspeed;
        if (Camera_fov<20.0) then Camera_fov=20.0;Camera_zoomin = false; end
        chiGraphicsCameraFOV(Camera_fov);
    end
    
    if (Camera_zoomout) then
        Camera_fov=Camera_fov+Camera_zoomspeed;
        if (Camera_fov>40.0) then Camera_fov=40.0; Camera_zoomout = false;end
        chiGraphicsCameraFOV(Camera_fov);
    end
    
    if ((Camera_left)                                      ) then chiGraphicsTranslateCamera(-Camera_movespeed,0,0); end
    if ((Camera_right)                                     ) then chiGraphicsTranslateCamera( Camera_movespeed,0,0); end
    if ((Camera_forward)  and (Camera_type=="Perspective") ) then chiGraphicsTranslateCamera(0, Camera_movespeed,0); end
    if ((Camera_backward) and (Camera_type=="Perspective") ) then chiGraphicsTranslateCamera(0,-Camera_movespeed,0); end
    if ((Camera_forward)  and (Camera_type=="Orthographic")) then Camera_orthoWidth=Camera_orthoWidth-Camera_movespeed; chiGraphicsCameraOrthoWidth(Camera_orthoWidth); end
    if ((Camera_backward) and (Camera_type=="Orthographic")) then Camera_orthoWidth=Camera_orthoWidth+Camera_movespeed; chiGraphicsCameraOrthoWidth(Camera_orthoWidth); end
    if ((Camera_up)                                        ) then chiGraphicsTranslateCamera(0, 0.0, Camera_movespeed); end
    if ((Camera_down)                                      ) then chiGraphicsTranslateCamera(0, 0.0,-Camera_movespeed); end

    Camera_x,Camera_y,Camera_z=chiGraphicsGetCameraPosition();
    
    if (Camera_x>Camera_xmax) then chiGraphicsTranslateCamera(-Camera_movespeed, 0.0, 0.0); end
    if (Camera_x<Camera_xmin) then chiGraphicsTranslateCamera( Camera_movespeed, 0.0, 0.0); end
    if (Camera_y>Camera_ymax) then chiGraphicsTranslateCamera( 0.0,-Camera_movespeed, 0.0); end
    if (Camera_y<Camera_ymin) then chiGraphicsTranslateCamera( 0.0, Camera_movespeed, 0.0); end
    if (Camera_z>Camera_zmax) then chiGraphicsTranslateCamera( 0.0, 0.0,-Camera_movespeed); end
    if (Camera_z<Camera_zmin) then chiGraphicsTranslateCamera( 0.0, 0.0, Camera_movespeed); end
    
    Camera_x,Camera_y,Camera_z=chiGraphicsGetCameraPosition();
end