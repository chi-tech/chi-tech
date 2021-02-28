mfdDisplayPanel=PanelClass.New("PreviewPanel");

mfdDisplayPanel.zDepth=1.1;
mfdDisplayPanel.xmin=1000;
mfdDisplayPanel.xmax=1200;
mfdDisplayPanel.ymin=1000;
mfdDisplayPanel.ymax=1200;
mfdDisplayPanel.fixedHeight=true;
mfdDisplayPanel.fixedWidth=true;
mfdDisplayPanel.AddBoundaryLock(mfdDisplayPanel,"NORTH",10,panels[7],true)
mfdDisplayPanel.AddBoundaryLock(mfdDisplayPanel,"WEST",10,panels[7],true)


dRegisterFeature(mfdDisplayPanel);