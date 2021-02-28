print("Creating newForm")
selectPhysicsModelsForm=FormClass.New("Select Physics Models",true)
dRegisterFeature(selectPhysicsModelsForm);
mtemp=selectPhysicsModelsForm

currentScene=chiGetScene();
chiBindScene(mtemp.sceneNumber);

mtemp.mainPanel = PanelClass.New("Main");

mtemp.RegisterFeature(mtemp,mtemp.mainPanel);
    
mtemp.mainPanel.AddBoundaryLock(mtemp.mainPanel,"NORTH");
mtemp.mainPanel.AddBoundaryLock(mtemp.mainPanel,"SOUTH");
mtemp.mainPanel.AddBoundaryLock(mtemp.mainPanel,"EAST");
mtemp.mainPanel.AddBoundaryLock(mtemp.mainPanel,"WEST");

--ambient=0.5;
--chiMaterialSetProperty(mtemp.mainPanel.matlNum,CHI_DIFFUSE_COLOR,ambient,ambient,ambient,1.0);

--newForm.Show(newForm)

--################################################ Labels
mtemp.labelStack={}
mtemp.checkBoxStack={}

--=============================================== Point Kinetics
lvl=1
mtemp.labelStack[lvl]=LabelClass.New("Reactor Point Kinetics");
mtemp.labelStack[lvl].SetProperty(mtemp.labelStack[lvl],"Master",mtemp.mainPanel);
mtemp.labelStack[lvl].SetProperty(mtemp.labelStack[lvl],"Float",false);
mtemp.labelStack[lvl].xSize=250;
mtemp.checkBoxStack[lvl] = CheckBoxClass.New(mtemp.labelStack[lvl].name.."_Checkbox")
mtemp.checkBoxStack[lvl].SetProperty(mtemp.checkBoxStack[lvl],"Master",mtemp.mainPanel);
mtemp.checkBoxStack[lvl].SetProperty(mtemp.checkBoxStack[lvl],"Float",false);
mtemp.checkBoxStack[lvl].xSize=20;
mtemp.RegisterFeature(mtemp,mtemp.labelStack[lvl])
mtemp.RegisterFeature(mtemp,mtemp.checkBoxStack[lvl]);

--=============================================== Smooth Particle Hydrodynamics
lvl=2
mtemp.labelStack[lvl]=LabelClass.New("Smooth Particle Hydrodynamics (SPH)");
mtemp.labelStack[lvl].SetProperty(mtemp.labelStack[lvl],"Master",mtemp.mainPanel);
mtemp.labelStack[lvl].SetProperty(mtemp.labelStack[lvl],"Float",false);
mtemp.labelStack[lvl].xSize=250;
mtemp.checkBoxStack[lvl] = CheckBoxClass.New(mtemp.labelStack[lvl].name.."_Checkbox")
mtemp.checkBoxStack[lvl].SetProperty(mtemp.checkBoxStack[lvl],"Master",mtemp.mainPanel);
mtemp.checkBoxStack[lvl].SetProperty(mtemp.checkBoxStack[lvl],"Float",false);
mtemp.checkBoxStack[lvl].xSize=20;
mtemp.RegisterFeature(mtemp,mtemp.labelStack[lvl])
mtemp.RegisterFeature(mtemp,mtemp.checkBoxStack[lvl]);


--############################################### Buttons
mtemp.applyButton = ButtonClass.New("Apply");
mtemp.RegisterFeature(mtemp,mtemp.applyButton);

func = function (this)
    --for k=1,1 do
        if (selectPhysicsModelsForm.checkBoxStack[1]~=nil) then
            if (selectPhysicsModelsForm.checkBoxStack[1].selected) then
                newPK=ReactorDynamics1DClass.New("Physics")
                chinTimeControl.AddSlave(chinTimeControl,newPK)
                print("PK Physics created")
            end
        end
        if (selectPhysicsModelsForm.checkBoxStack[2]~=nil) then
            if (selectPhysicsModelsForm.checkBoxStack[2].selected) then
                if (SPHCount==0) then   --Limit creation to 1
                    newSPH=SPHClass.New("SPH");
                    chinTimeControl.AddSlave(chinTimeControl,newSPH)
                end
            end
        end
    --end
end
mtemp.applyButton.CustomButtonDown=func



chiBindScene(currentScene);