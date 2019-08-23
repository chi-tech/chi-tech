Camera_alpha        = 0;        -- camera orientation on the x axis
Camera_azimuth      = -80.0;    -- camera orientation on the y axis
Camera_tsi          = 0;        -- something
Camera_x            = 0.0;      -- x pos. 
Camera_y            = -5.0;    -- y pos.
Camera_z            =  2.0;      -- z pos.
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
Camera_refPoint     = {0,0,0};


chiGraphicsOrientCamera(0,Camera_alpha,Camera_azimuth,Camera_tsi);
chiGraphicsPositionCamera(Camera_x,Camera_y,Camera_z);
spec=1.0;
xspec=0.5;
yspec=0.0;

function CameraControl()
    
    if (WM_CTLBUTTONDN.occured) then Camera_lock1open = true; end
    if (WM_MBUTTONDOWN.occured) then Camera_lock2open = true; end
    if (WM_RBUTTONDOWN.occured) then Camera_lock1open = true; end
    
    if (WM_CTLBUTTONUP.occured) then Camera_lock1open = false; end
    if (WM_MBUTTONUP.occured) then Camera_lock2open = false; end
    if (WM_RBUTTONUP.occured) then Camera_lock1open = false; end

    
    
    --================================= View panning
    if ((Camera_lock1open) and (Camera_lock2open) and (WM_MOUSEMOVE.occured)) then
        xr,yr,zr = chi3DPointToScreen(Camera_refPoint[1],Camera_refPoint[2],Camera_refPoint[3]);
        wx1,wy1,wz1 = chiScreenTo3DPoint(WM_MOUSEMOVE.iPar0,WM_MOUSEMOVE.iPar1,zr);
        wx2,wy2,wz2 = chiScreenTo3DPoint(WM_MOUSEMOVE.iPar2,WM_MOUSEMOVE.iPar3,zr);
        
        mouse_dx = WM_MOUSEMOVE.iPar0 - WM_MOUSEMOVE.iPar2;
        mouse_dy = WM_MOUSEMOVE.iPar1 - WM_MOUSEMOVE.iPar3;
        
        ratio=math.abs(mouse_dx)/(math.abs(mouse_dx)+math.abs(mouse_dy));

        sum1=     (wx1-wx2)*(wx1-wx2);
        sum1=sum1+(wy1-wy2)*(wy1-wy2);
        sum1=sum1+(wz1-wz2)*(wz1-wz2);
        sum2=math.sqrt(sum1);
        translation_x=sum2*2*ratio;
        translation_z=sum2*2*(1-ratio);

        if     (mouse_dx>0) then
            chiGraphicsTranslateCamera(-translation_x,0,0);
        elseif (mouse_dx<0) then 
            chiGraphicsTranslateCamera(translation_x,0,0);
        end
        if     (mouse_dy>0) then
            chiGraphicsTranslateCamera(0,0, translation_z);
        elseif (mouse_dy<0) then        
            chiGraphicsTranslateCamera(0,0,-translation_z);
        end
    end
    
    --================================= View rotation
    if ((not Camera_lock1open) and (Camera_lock2open) and (WM_MOUSEMOVE.occured)) then
        mouse_dx = WM_MOUSEMOVE.iPar0 - WM_MOUSEMOVE.iPar2;
        mouse_dy = WM_MOUSEMOVE.iPar1 - WM_MOUSEMOVE.iPar3;
        Camera_alpha=Camera_alpha-0.2*mouse_dx;
        Camera_azimuth=Camera_azimuth-0.2*mouse_dy;
        chiGraphicsOrientCamera(0,Camera_alpha,Camera_azimuth,Camera_tsi);
    end
    
    --================================= View zoom
    if (WM_MOUSEWHEEL.occured) then
        mouse_dz = -WM_MOUSEWHEEL.iPar0*Camera_zoomspeed/120.0;
        if     (mouse_dz>0) then
            chiGraphicsTranslateCamera(0, Camera_movespeed,0);
        elseif (mouse_dz<0) then                          
            chiGraphicsTranslateCamera(0,-Camera_movespeed,0);
        end
    end
    
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