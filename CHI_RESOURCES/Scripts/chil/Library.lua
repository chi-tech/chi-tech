chiLibraryRoot = "CHI_RESOURCES/Scripts/chil/"
print("Checkpoint1")
--===================================== UtilityFunctions
dofile(chiLibraryRoot.."Callbacks.lua")
dofile(chiLibraryRoot.."../chi_Nteractor/chin_GlobalFunctions/chin_globalfs_getscriptpath.lua")
dofile(chiLibraryRoot.."SelectionStack.lua")


--===================================== Cameras
dofile(chiLibraryRoot.."Cameras/cameras.lua")
dofile(chiLibraryRoot.."Cameras/CHI_CAMERA_THIRDPERSON.lua")
dofile(chiLibraryRoot.."Cameras/CHI_CAMERA_ORTHOWINDOW.lua")
dofile(chiLibraryRoot.."Cameras/CHI_CAMERA_ORTHODISPLAYER.lua")
dofile(chiLibraryRoot.."Cameras/CHI_CAMERA_REVOLVER.lua")

--===================================== Default Material
chiDefaultMat  = chiMaterialCreate()


--===================================== Default light


--===================================== World items
CHI_WORLD_DEFAULTFLOOR = chiLibraryRoot.."WorldItems/CHI_WORLD_DEFAULTFLOOR.lua"

--===================================== Forms
dofile(chiLibraryRoot.."FormElements/forms.lua")
dofile(chiLibraryRoot.."FormElements/forms_01_vsplit.lua")
dofile(chiLibraryRoot.."FormElements/forms_03_panel.lua")



--===================================== Defining the main function
function main()
    chilCallbacks.Execute()
end

print("Library loaded")