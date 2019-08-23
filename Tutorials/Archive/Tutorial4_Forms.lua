--===================================== Load library
dofile(CHI_LIBRARY)

chiSetWindowProperties(400,400,200,200)

baseScene,displayer2D = chiGetScene();
camera2D = chilCreateOrthoWindowCamera("camera2D")
chilCameraOrganizer.AddCamera(camera2D);

chilForms_CreateVSplit()

newForm = chilForms_CreateExternalForm("New Form");

