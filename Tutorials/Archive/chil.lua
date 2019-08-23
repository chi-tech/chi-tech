--===================================== Load library
if (CHI_LIBRARY~=nil) then
  dofile(CHI_LIBRARY)
else
  print("ERROR: CHI_LIBRARY not found")
  return;
end

mainCamera = chilCreateOrthoWindowCamera("MainCamera")
chilCameraOrganizer.AddCamera(mainCamera);

chilForms_CreatePanel("TopPanel")
chilForms_CreateVSplit("Middl")