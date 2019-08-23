

--=============================================== Menubar
menuItem={}
menuItem[1]=MenuClass.New("File");
menuItem[2]=MenuClass.New("Edit");
menuItem[3]=MenuClass.New("Tools");

menuItem[1].paddingLeft=5;

for k=1,3 do
    menuItem[k].SetProperty(menuItem[k],"Master",panels[1]);
    menuItem[k].float=false
    menuItem[k].xSize=50;
    menuItem[k].ySize=28;
    menuItem[k].panel1.cutWindow=1;
    menuItem[k].panel1.cutPanel=panels[4];
    dRegisterFeature(menuItem[1]);
end

size    = 20;
padLeft = 20;
padBot  = 4;

menuItem[1].subItems={}

newLabel=menuItem[1].AddListItem(menuItem[1],"New Simulation");
newLabel.paddingLeft=padLeft;
newLabel.paddingBot=padBot;
newLabel.ySize=size+10;
dRegisterFeature(newLabel)
menuItem[1].subItems[1]=newLabel;


newLabel=menuItem[1].AddListItem(menuItem[1],"Open");
newLabel.paddingLeft=padLeft;
newLabel.paddingBot=padBot;
newLabel.ySize=size;
dRegisterFeature(newLabel)
menuItem[1].subItems[2]=newLabel;

newLabel=menuItem[1].AddListItem(menuItem[1],"Save");
newLabel.paddingLeft=padLeft;
newLabel.paddingBot=padBot;
newLabel.ySize=size;
dRegisterFeature(newLabel)
menuItem[1].subItems[3]=newLabel;

newLabel=menuItem[1].AddListItem(menuItem[1],"Save As");
newLabel.paddingLeft=padLeft;
newLabel.paddingBot=padBot;
newLabel.ySize=size;
dRegisterFeature(newLabel)
menuItem[1].subItems[4]=newLabel;

newLabel=menuItem[1].AddListItem(menuItem[1],"Import Object");
newLabel.paddingLeft=padLeft;
newLabel.paddingBot=padBot;
newLabel.ySize=size;
dRegisterFeature(newLabel)
menuItem[1].subItems[5]=newLabel;

--newLabel=menuItem[1].AddListItem(menuItem[1],"New Part");
--newLabel.paddingLeft=padLeft;
--newLabel.paddingBot=padBot;
--newLabel.ySize=size;
--dRegisterFeature(newLabel)
folderPath = chinGetScriptPath()
dofile(folderPath.."/Menu/MainMenu_File_00.lua")

