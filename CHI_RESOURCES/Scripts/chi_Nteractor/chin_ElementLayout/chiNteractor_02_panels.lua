--=============================================== Panels
panels                      = {};
panels[0]                   = PanelClass.New("Background");
panels[1]                   = PanelClass.New("Menu Bar");
panels[2]                   = PanelClass.New("Tool Bar");
panels[3]                   = PanelClass.New("Tree View");
panels[4]                   = PanelClass.New("Scene View");
panels[5]                   = PanelClass.New("Properties");
panels[6]                   = PanelClass.New("Console");
panels[7]                   = PanelClass.New("Information Pane");
panels[8]                   = PanelClass.New("Status Bar");

dRegisterFeature(panels[0]);
dRegisterFeature(panels[1]);
dRegisterFeature(panels[2]);
dRegisterFeature(panels[3]);
dRegisterFeature(panels[4]);
dRegisterFeature(panels[5]);
dRegisterFeature(panels[6]);
dRegisterFeature(panels[7]);
dRegisterFeature(panels[8]);

ambient=0.97;
chiMaterialSetProperty(panels[1].matlNum,CHI_DIFFUSE_COLOR,ambient,ambient,ambient,1.0);
chiMaterialSetProperty(panels[2].matlNum,CHI_DIFFUSE_COLOR,ambient,ambient,ambient,1.0);

--chiMaterialSetProperty(panels[7].matlNum,CHI_DIFFUSE_COLOR,0.25,0.25,1.0,1.0);

--

panels[0].AddBoundaryLock(panels[0],"NORTH");
panels[0].AddBoundaryLock(panels[0],"SOUTH");
panels[0].AddBoundaryLock(panels[0],"EAST");
panels[0].AddBoundaryLock(panels[0],"WEST");
chiMaterialSetProperty(panels[0].matlNum,"Ambient",0.5,0.5,0.5,1.0);
panels[0].zDepth=0.99;

--
panels[1].xmin              = 0;
panels[1].xmax              = chinGlobal.dwindowxsize;
panels[1].ymin              = chinGlobal.dwindowysize-30
panels[1].ymax              = chinGlobal.dwindowysize;
panels[1].AddBoundaryLock(panels[1],"NORTH");
panels[1].AddBoundaryLock(panels[1],"EAST");
panels[1].fixedHeight       = true;
panels[1].Redraw(panels[1]);

--
panels[2].xmin              = 0;
panels[2].xmax              = chinGlobal.dwindowxsize;
panels[2].ymin              = panels[1].ymin-30;
panels[2].ymax              = panels[1].ymin;
panels[2].AddBoundaryLock(panels[2],"NORTH",0,panels[1]);
panels[2].AddBoundaryLock(panels[2],"EAST");
panels[2].fixedHeight       = true;
panels[2].Redraw(panels[2]);

--              
panels[3].xmin              = 0;
panels[3].xmax              = 200;
panels[3].ymin              = 200;
panels[3].ymax              = panels[2].ymin-2;
panels[3].AddBoundaryLock(panels[3],"NORTH",0,panels[2]);
panels[3].AddScrollbar(panels[3],10.0,0,true);
panels[3].Redraw(panels[3]);

--              
panels[4].xmin              = 202;
panels[4].xmax              = chinGlobal.dwindowxsize;
panels[4].ymin              = 200;
panels[4].ymax              = panels[2].ymin-2;
panels[4].AddBoundaryLock(panels[4],"NORTH",0,panels[2]);
panels[4].AddBoundaryLock(panels[4],"EAST");
panels[4].Redraw(panels[4]);
--
panels[4].AttachWindow(panels[4],1);

--
panels[5].xmin              = 0;
panels[5].xmax              = 200;
panels[5].ymin              = 20;
panels[5].ymax              = 198;
panels[5].AddScrollbar(panels[5],10.0,0,true);
panels[5].Redraw(panels[5]);
--panels[5].AddBoundaryLock(panels[5],"SOUTH",panels[7]);

--
panels[6].xmin              = 204 ;
panels[6].xmax              = chinGlobal.dwindowxsize;
panels[6].ymin              = 20;
panels[6].ymax              = 198;
panels[6].AddScrollbar(panels[6],10.0,0,true);
--panels[6].widthMin=700;
--panels[6].AddBoundaryLock(panels[6],"EAST");
--panels[6].AddBoundaryLock(panels[6],"SOUTH",panels[7]);
panels[6].Redraw(panels[6]);

--
panels[7].ymin              = 20;
panels[7].xmin              = 960;
panels[7].AddBoundaryLock(panels[7],"EAST");

--
panels[8].ymax              = 20;
panels[8].AddBoundaryLock(panels[8],"SOUTH");
panels[8].AddBoundaryLock(panels[8],"EAST");
panels[8].AddBoundaryLock(panels[8],"WEST");

