--#############################--
--         Tutorial 1
--  Basic window, text and line creation
--#############################--

chiSetWindowProperties(400,400,200,200)
line3Dindex = chi3DLineCreate("3DLine")
zheight=0.5;
chi3DLineChangeColor(line3Dindex,0,0,1,1);      --Blue line
chi3DLineAddVertex(line3Dindex,  0,  0, zheight);
chi3DLineAddVertex(line3Dindex,  2,  0, zheight);
chi3DLineAddVertex(line3Dindex,  1,  2, zheight);
chi3DLineAddVertex(line3Dindex,  0,  0, zheight);





chiWindowCreate("Bob",false)
chiBindScene(1)
chiSetWindowProperties(800,400,800,200)
line3Dindex2 = chi3DLineCreate("3DLine2")
chi3DLineChangeColor(line3Dindex2,1,0,0,1);      --Red line
chi3DLineAddVertex(line3Dindex2,  0,  0, 0);
chi3DLineAddVertex(line3Dindex2,  3,  0, 0);
chi3DLineAddVertex(line3Dindex2,  1,  2, 0);
chi3DLineAddVertex(line3Dindex2,  0,  0, 0);
deferWindowSetupDone=true;

dofile(CHI_WORLD_DEFAULTFLOOR)

mainCamera = chiCreateThirdPersonCamera("MainCamera")
chiSetActiveCamera(mainCamera)






