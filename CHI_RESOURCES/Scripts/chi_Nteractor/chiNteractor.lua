chinBaseDir = "CHI_RESOURCES/Scripts/chi_Nteractor/";

dofile(chinBaseDir.."chiNteractor_00_global.lua")
dofile(chinBaseDir.."chiNteractor_02_icons.lua")
dofile(chinBaseDir.."chin_Form/chin_form.lua")

if (CUSTOM_LAYOUT==nil) then
    --chiSetSceneUpdateMode(1);
    dofile(chinBaseDir.."/chin_ElementLayout/chiNteractor_00_layout.lua")
    chiSetWindowProperties("MAXIMIZED")
end

chiBindScene(1)
if (CUSTOM_LAYOUT==nil) then
    --chiSetSceneUpdateMode(1);
    dofile(chinBaseDir.."chin_3dscene/chin_3dscene.lua")
end







function d_main()

    chiBindScene(0);
    chiRequestSceneRefresh();
    
    d_events(); --Process events

    if (not (CustomProcessing==nil)) then 
        CustomProcessing(); 
    end
    
    d_filter(); --Filter events to not go to subscenes

    chiBindScene(1)
end




