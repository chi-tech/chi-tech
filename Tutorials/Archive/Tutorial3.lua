--===================================== Load library
dofile(CHI_LIBRARY)

chiSetWindowProperties(400,400,200,200)

baseScene,displayer2D = chiGetScene();
camera2D = chilCreateOrthoWindowCamera("camera2D")
chilCameraOrganizer.AddCamera(camera2D);

--text1 = chiTextCreate("Text");
--chiTextSetProperty(text1,1,"Hello World");
--chiTextSetProperty(text1,2,200,200,0.0);



displayer3D = chiDisplayerCreate();
chiBindScene(baseScene,displayer3D);

dofile(CHI_WORLD_DEFAULTFLOOR);

camera3D = chilCreateThirdPersonCamera("camera3D")
chilCameraOrganizer.AddCamera(camera3D);
--
--line3Dindex = chi3DLineCreate("3DLine")
--chi3DLineChangeColor(line3Dindex,0,0,1,1);      --Blue line
--chi3DLineAddVertex(line3Dindex,  0,  0, 0);
--chi3DLineAddVertex(line3Dindex,  2,  0, 0);
--chi3DLineAddVertex(line3Dindex,  1,  2, 0);
--chi3DLineAddVertex(line3Dindex,  0,  0, 0);
--
--
--surf=chiObjectLoadSurface("CHI_RESOURCES/TestObjects/Monkey.obj");
--monkeyIndex=chiObjectCreate("Monkey");
--chiObjectAddSurface(monkeyIndex,surf);
--
--monkeyTrans = chiTransformCreate()
----chiTransformSetRotationPoint(monkeyTrans,-100.0,-100.0,0.0)
----chiTransformSetRotation(monkeyTrans,0.0,0.0,15.0);
----chiTransformSetTranslation(monkeyTrans,100.0,100.0,0.0);
----chiTransformSetScale(monkeyTrans,50.0,50.0,50.0);
--
--chiObjectSetProperty(monkeyIndex,"Transform",monkeyTrans);





