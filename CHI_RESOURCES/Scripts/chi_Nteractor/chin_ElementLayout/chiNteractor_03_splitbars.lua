--=============================================== Split bars
splitbars={}

splitbars[1]=SplitBarClass.New("SplitV1")
splitbars[2]=SplitBarClass.New("SplitH1")
splitbars[3]=SplitBarClass.New("SplitH2")
splitbars[4]=SplitBarClass.New("SplitV2")

dRegisterFeature(splitbars[1]);
dRegisterFeature(splitbars[2]);
dRegisterFeature(splitbars[3]);
dRegisterFeature(splitbars[4]);

splitbars[1].gripWidth=8;
splitbars[2].gripWidth=8;
splitbars[3].gripWidth=8;
splitbars[4].gripWidth=8;

splitbars[1].ymin=0; splitbars[1].xmin=250;

splitbars[1].AddSlave(splitbars[1],panels[3],"EAST",0)
splitbars[1].AddSlave(splitbars[1],panels[4],"WEST", 4)
splitbars[1].AddSlave(splitbars[1],panels[5],"EAST",0)
splitbars[1].AddSlave(splitbars[1],panels[6],"WEST", 4)

splitbars[2].ymin=250; splitbars[2].xmin=0; 
splitbars[2].vertical=false;

splitbars[2].AddSlave(splitbars[2],panels[3],"SOUTH",4)
splitbars[2].AddSlave(splitbars[2],panels[5],"NORTH",0)

splitbars[3].ymin=250; splitbars[3].xmin=0;
splitbars[3].vertical=false;

splitbars[3].AddSlave(splitbars[3],panels[4],"SOUTH",4)
splitbars[3].AddSlave(splitbars[3],panels[6],"NORTH",0)
splitbars[3].AddSlave(splitbars[3],panels[7],"NORTH",0)

splitbars[4].xmin=950;
splitbars[4].AddSlave(splitbars[4],panels[6],"EAST",0)
splitbars[4].AddSlave(splitbars[4],panels[7],"WEST",4)

splitbars[1].Redraw(splitbars[1]);
splitbars[2].Redraw(splitbars[2]);
splitbars[3].Redraw(splitbars[3]);
splitbars[4].Redraw(splitbars[4]);