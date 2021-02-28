


function chilForms_CreatePanel(name)
  local newPanel = {}
  newPanel.name                = name;
  newPanel.sceneScope          = -1;
  newPanel.displayerScope      = -1;
  --.camera
  newPanel.callbackReference   = newPanel;
  newPanel.callBackObject      = nil;

  --newPanel.initialized         = false;

  newPanel.anchors = {left=nil,right=nil,top=nil,bottom=nil }
  newPanel.parent  = nil;
  newPanel.children = {itemCount=0 }
  newPanel.selected = false;

  newPanel.position = {x=100,y=100};
  newPanel.size     = {width=200,height=200};

  newPanel.lastObject          = -1;

  newPanel.sceneScope, newPanel.displayerScope = chiGetScene();

  newPanel.divider = chi3DLineCreate("");

  chi3DLineAddVertex(newPanel.divider,newPanel.position.x,newPanel.position.y, 0.0);
  chi3DLineAddVertex(newPanel.divider,newPanel.position.x+newPanel.size.width,newPanel.position.y, 0.0);
  chi3DLineAddVertex(newPanel.divider,newPanel.position.x+newPanel.size.width,newPanel.position.y+newPanel.size.height, 0.0);
  chi3DLineAddVertex(newPanel.divider,newPanel.position.x,newPanel.position.y+newPanel.size.height, 0.0);
  chi3DLineAddVertex(newPanel.divider,newPanel.position.x,newPanel.position.y, 0.0);
  chi3DLineSetStipple(newPanel.divider,false,1,5.0);

  --====================================================== Callback function
  newPanel.callbackFunction = function (this)
    --=================================== Mouse move event

  end

  --====================================================== Callback object
  newPanel.callBackObject = chilCallbacks.MakeCallback(newPanel);
  chilCallbacks.PushCallback(newPanel.callBackObject);

  --====================================================== Select method
  newPanel.Select = function (this)
    this.selected = true;
  end

  --====================================================== DeSelect method
  newPanel.DeSelect = function (this)
    this.selected = false;
  end

  return newPanel;
end