--######################################################### Add slave
function ManipulatorClass.AddObjectSlave(this,theslave)
    this.objectSlaveCount = this.objectSlaveCount + 1;
    k=this.objectSlaveCount;
    
    this.objectSlaves[k]=theslave;
    
    x,y,z = chiTransformGet(theslave.obj1Tra);
    chiTransformSetTranslation(this.manipTransform,x,y,z);
end

--######################################################### Add slave
function ManipulatorClass.RemoveObjectSlave(this,theslave)
    
    
end

--######################################################### Clear slaves
function ManipulatorClass.RemoveAllObjectSlaves(this)
    this.objectSlaveCount = 0;
end

--########################################################## Master translate
function ManipulatorClass.MasterTranslate(this,dx,dy,dz)
    for k=1,this.objectSlaveCount do
        theSlave=this.objectSlaves[k];
        chiTransformSetTranslation(theSlave.obj1Tra,dx,dy,dz);
    end
    chiTransformSetTranslation(this.manipTransform,dx,dy,dz);
end
