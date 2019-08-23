--######################################################### Constr
PARTICLE_POSITION=1
PARTICLE_VELOCITY=2
PARTICLE_DENSITY =3
PARTICLE_SMOOTHINGLENGTH=4
PARTICLE_DPDRHO=5
PARTICLE_PRESSURE=6
PARTICLE_MASS=7
PARTICLE_FIXED=8
function SPHClass.New(name)
    local this = setmetatable({},SPHClass);
    SPHCount = SPHCount + 1;
    
    local objName = name.."_Domain";
    this.objName = objName;
    
    this.sphObj = chiSPHCreate();
    currentScene=chiGetScene();
    chiBindScene(1);
    this.domainSurf = chiLoadSurface("CHI_RESOURCES/TestObjects/CubeDomain.obj");
    this.domainObj = chiObjectCreate(objName); 
    chiObjectAddSurface(objName,this.domainSurf);
    
    this.matlNum = chiMaterialCreate(objName .. "_Material");
    ambient = 1.0;
    chiMaterialSetProperty(this.matlNum,"Diffuse",ambient,ambient,ambient,1.0);
    chiObjectSetProperty(objName,"Material",objName .. "_Material");
    chiObjectSetProperty(objName,"Wireframe",true);
    
    --================ Create N particles
    this.N = 10;
    dmin=-1;
    dmax= 1;
    deltaD = (dmax-dmin)/(this.N+1);
    psize = 15;
    chi3DPointCreate(objName.."_Points");
    
    counter=-1;
    for i=1,this.N do
        for j=1,this.N do
            for k=1,this.N do
                counter=counter+1;
                x=deltaD*i+dmin;
                y=deltaD*j+dmin;
                z=deltaD*k+1.0+dmin;
                --print(PARTICLE_POSITION,counter,i,j,k)
                chiSPHSetProperty(this.sphObj,PARTICLE_POSITION,counter,x,y,z);
                chi3DPointAddVertex(objName.."_Points",x,y,z);
                chi3DPointChangeVertex(objName.."_Points",counter,x,y,z,0,0,0,1,psize);
                
            end
        end
    end
    
    print("SPH Created as object: ",objName)
    chiBindScene(currentScene);
    return this;
end

function SPHClass.Initialize(this)
    --================ Create N particles
    this.N = 10;
    dmin=-1;
    dmax= 1;
    deltaD = (dmax-dmin)/(this.N+1);
    psize = 15;
    print("Initializing SPH")
    counter=-1;
    
    for i=1,this.N do
        for j=1,this.N do
            for k=1,this.N do
                counter=counter+1;
                x=deltaD*i+dmin;
                y=deltaD*j+dmin;
                z=deltaD*k+1.0+dmin;
                
                
                chi3DPointChangeVertex(this.objName.."_Points",counter,x,y,z,0,0,0,1,psize);
                chiSPHSetProperty(this.sphObj,PARTICLE_POSITION,counter,x,y,z);
                chiSPHSetProperty(this.sphObj,PARTICLE_VELOCITY,counter,0,0,0);
                chiSPHSetProperty(this.sphObj,PARTICLE_DENSITY,counter,1000.0);
                chiSPHSetProperty(this.sphObj,PARTICLE_DPDRHO,counter,2000.0);
                chiSPHSetProperty(this.sphObj,PARTICLE_PRESSURE,counter,0.0);
                chiSPHSetProperty(this.sphObj,PARTICLE_MASS,counter,1000.0*deltaD*deltaD*deltaD);
                chiSPHSetProperty(this.sphObj,PARTICLE_SMOOTHINGLENGTH,counter,deltaD);
                if (k==1) then
                    chiSPHSetProperty(this.sphObj,PARTICLE_FIXED,counter,true);
                end
                
            end
        end
    end
    chiSPHInitialize();
    print("SPH Initialized")
end

function SPHClass.Step(this,dt)
    chiSPHStep(this.sphObj,dt);
    counter=-1;
    for i=1,this.N do
        for j=1,this.N do
            for k=1,this.N do
                counter=counter+1;
                x,y,z=chiSPHQuery(this.sphObj,PARTICLE_POSITION,counter);
                chi3DPointChangeVertex(this.objName.."_Points",counter,x,y,z,0,0,0,1,psize);
            end
        end
    end
end