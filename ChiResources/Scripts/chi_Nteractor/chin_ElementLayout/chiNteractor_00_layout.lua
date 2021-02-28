dofile(chinBaseDir.."/chin_ElementLayout/chiNteractor_00_updatextechelements.lua")
dofile(chinBaseDir.."/chin_ElementLayout/chiNteractor_01_windowsetup.lua")
dofile(chinBaseDir.."/chin_ElementLayout/chiNteractor_02_panels.lua")
dofile(chinBaseDir.."/chin_ElementLayout/chiNteractor_03_splitbars.lua")
dofile(chinBaseDir.."/chin_ElementLayout/chiNteractor_04_treeview.lua")
dofile(chinBaseDir.."/chin_ElementLayout/chiNteractor_05_console.lua")
dofile(chinBaseDir.."/chin_ElementLayout/chiNteractor_06_gridview.lua")
dofile(chinBaseDir.."/chin_ElementLayout/chiNteractor_07_menubar.lua")
dofile(chinBaseDir.."/chin_ElementLayout/chiNteractor_08_mfd.lua")
dofile(chinBaseDir.."/chin_ElementLayout/chiNteractor_09_contextmenus.lua")
dofile(chinBaseDir.."/chin_ElementLayout/chiNteractor_10_toolbar.lua")
dofile(chinBaseDir.."/chin_ElementLayout/chiNteractor_11_subforms.lua")






chiBindScene(0);

--######################################################### Custom processing
function CustomProcessing()

    UpdateXTechElements();
    
    TreeViewUpdate();
    

end


