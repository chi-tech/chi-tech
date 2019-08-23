--#############################--
--         Tutorial 1
--  Basic window, text and line creation
--#############################--
print("Tutorial 1")
--===================================== Load library
if (CHI_LIBRARY~=nil) then
  dofile(CHI_LIBRARY)
else
  print("ERROR: CHI_LIBRARY not found")
  return;
end

print("Tutorial 1 checkpoint")
baseScene,displayer0 = chiGetScene();
mainCamera = chilCreateOrthoWindowCamera("MainCamera")
chilCameraOrganizer.AddCamera(mainCamera);

--===================================== Window setup
--Setting window properties xpos,ypos,xwidth,ywidth
chiSetWindowProperties(400,400,200,200)

--===================================== Text creation
text1 = chiTextCreate("Text");
chiTextSetProperty(text1,1,"Hello World");
chiTextSetProperty(text1,2,100,100,0);

--===================================== Line creation in world space
--The default camera is an orthographic camera looking down the z-axis and with
--its up-vector on the y-axis. It is located +10 units off the z=0 plane and
--has a frustum of 0.1 to 999 units. Its x- and y-units are default x=[-5,5] and
--y=x/aspectratio and can be changed using chiGraphicsCameraOrthoWidth(10)
line3Dindex = chi3DLineCreate("3DLine")
chi3DLineChangeColor(line3Dindex,0,0,1,1);      --Blue line
chi3DLineAddVertex(line3Dindex,  0,  1, 0);
chi3DLineAddVertex(line3Dindex,  200,  1, 0);
chi3DLineAddVertex(line3Dindex,  100,  200, 0);
chi3DLineAddVertex(line3Dindex,  0,  1, 0);


--displayer1 = chiDisplayerCreate();
--chiBindScene(baseScene,displayer1);
--
--chiDisplayerSetViewport(200,200,400,400);
--
--displayer1Camera = chilCreateOrthoDisplayerCamera("Displayer1Camera");
--chilCameraOrganizer.AddCamera(displayer1Camera);
--
--line3Dindex2 = chi3DLineCreate("3DLine2")
--chi3DLineChangeColor(line3Dindex2,1,0,0,1);      --Red line
--chi3DLineAddVertex(  line3Dindex2,  20,  1, 0);
--chi3DLineAddVertex(  line3Dindex2,  200,  1, 0);
--chi3DLineAddVertex(  line3Dindex2,  100,  200, 0);
--chi3DLineAddVertex(  line3Dindex2,  20,  1, 0);
