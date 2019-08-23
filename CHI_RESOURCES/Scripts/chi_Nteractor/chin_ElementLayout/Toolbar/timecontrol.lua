textNum = chiLoadTexture("CHI_RESOURCES/Textures/chin_Icons.png")
timeControl={}
timeControl.Button = {}
butNum=1
timeControl.Button[butNum] = ButtonClass.New("Initialize")
timeControl.Button[butNum].SetProperty(timeControl.Button[butNum],"Float",false)
timeControl.Button[butNum].SetProperty(timeControl.Button[butNum],"Master",panels[2])
timeControl.Button[butNum].SetProperty(timeControl.Button[butNum],"Text","")
timeControl.Button[butNum].xSize = 30;
--========================================================= Initialize button
timeControl.Button[butNum].ySize = 30;
chiMaterialSetProperty(timeControl.Button[butNum].matlNum,8,textNum);
chiTransformSetScale(timeControl.Button[butNum].obj1TTra,chinGlobal.Icons.iconScale,chinGlobal.Icons.iconScale,1.0);
dx = chinGlobal.Icons.iconShift*chinIconIni[1]
dy = chinGlobal.Icons.iconShift*chinIconIni[2]

chiTransformSetTranslation(timeControl.Button[butNum].obj1TTra, chinGlobal.Icons.iconTranslation +dx, chinGlobal.Icons.iconTranslation +dy, 0.0);
item=timeControl.Button[butNum]
function item.CustomButtonDown(this)
    chinTimeControl.Initialize(chinTimeControl)
end    

dRegisterFeature(timeControl.Button[butNum])

--========================================================= Step button
butNum=2
timeControl.Button[butNum] = ButtonClass.New("Step")
timeControl.Button[butNum].SetProperty(timeControl.Button[butNum],"Float",false)
timeControl.Button[butNum].SetProperty(timeControl.Button[butNum],"Master",panels[2])
timeControl.Button[butNum].SetProperty(timeControl.Button[butNum],"Text","")
timeControl.Button[butNum].xSize = 30;
timeControl.Button[butNum].ySize = 30;
chiMaterialSetProperty(timeControl.Button[butNum].matlNum,8,textNum);
chiTransformSetScale(timeControl.Button[butNum].obj1TTra,chinGlobal.Icons.iconScale,chinGlobal.Icons.iconScale,1.0);
dx = chinGlobal.Icons.iconShift*chinIconStep[1]
dy = chinGlobal.Icons.iconShift*chinIconStep[2]

chiTransformSetTranslation(timeControl.Button[butNum].obj1TTra, chinGlobal.Icons.iconTranslation +dx, chinGlobal.Icons.iconTranslation +dy, 0.0);
item=timeControl.Button[butNum]
function item.CustomButtonDown(this)
    chinTimeControl.Step(chinTimeControl)
end      

dRegisterFeature(timeControl.Button[butNum])

--========================================================= Run button
butNum=3
timeControl.Button[butNum] = ButtonClass.New("Run")
timeControl.Button[butNum].SetProperty(timeControl.Button[butNum],"Float",false)
timeControl.Button[butNum].SetProperty(timeControl.Button[butNum],"Master",panels[2])
timeControl.Button[butNum].SetProperty(timeControl.Button[butNum],"Text","")
timeControl.Button[butNum].xSize = 30;
timeControl.Button[butNum].ySize = 30;
chiMaterialSetProperty(timeControl.Button[butNum].matlNum,8,textNum);
chiTransformSetScale(timeControl.Button[butNum].obj1TTra,chinGlobal.Icons.iconScale,chinGlobal.Icons.iconScale,1.0);
dx = chinGlobal.Icons.iconShift*chinIconRun[1]
dy = chinGlobal.Icons.iconShift*chinIconRun[2]

chiTransformSetTranslation(timeControl.Button[butNum].obj1TTra, chinGlobal.Icons.iconTranslation +dx, chinGlobal.Icons.iconTranslation +dy, 0.0);
item=timeControl.Button[butNum]
function item.CustomButtonDown(this)
    chinTimeControl.running=true
end      

dRegisterFeature(timeControl.Button[butNum])

--========================================================= Stop button
butNum=4
timeControl.Button[butNum] = ButtonClass.New("Stop")
timeControl.Button[butNum].SetProperty(timeControl.Button[butNum],"Float",false)
timeControl.Button[butNum].SetProperty(timeControl.Button[butNum],"Master",panels[2])
timeControl.Button[butNum].SetProperty(timeControl.Button[butNum],"Text","")
timeControl.Button[butNum].xSize = 30;
timeControl.Button[butNum].ySize = 30;
chiMaterialSetProperty(timeControl.Button[butNum].matlNum,8,textNum);
chiTransformSetScale(timeControl.Button[butNum].obj1TTra,chinGlobal.Icons.iconScale,chinGlobal.Icons.iconScale,1.0);
dx = chinGlobal.Icons.iconShift*chinIconStop[1]
dy = chinGlobal.Icons.iconShift*chinIconStop[2]

chiTransformSetTranslation(timeControl.Button[butNum].obj1TTra, chinGlobal.Icons.iconTranslation +dx, chinGlobal.Icons.iconTranslation +dy, 0.0);
item=timeControl.Button[butNum]
function item.CustomButtonDown(this)
    chinTimeControl.running=false
end      

dRegisterFeature(timeControl.Button[butNum])


