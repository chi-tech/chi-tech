function ReactorDynamics1DClass.New(name)
    local this = setmetatable({},ReactorDynamics1DClass);
    ReactorDynamics1DCount = ReactorDynamics1DCount + 1;
    
    local objName="ReactorDynamics1D"..string.format("%03d",ReactorDynamics1DCount);
    this.name = name;
    
    --================================= Create the module
    this.moduleIndex = chiReactorDynamics1DCreate(name);
    
    
    
    
    
    return this;
end

function ReactorDynamics1DClass.Initialize(this)
    
    print("Initializing Point Kinetics")
    chiReactorDynamics1DInitialize(this.moduleIndex)
    chiReactorDynamics1DQuery(4,this.moduleIndex);
    print("Power",chinRDynamic1D[this.moduleIndex].initPower)
    print("Rho",chinRDynamic1D[this.moduleIndex].initRho)
    print("c1",chinRDynamic1D[this.moduleIndex].Precursor[1].ci)
    print("c2",chinRDynamic1D[this.moduleIndex].Precursor[2].ci)
    print("c3",chinRDynamic1D[this.moduleIndex].Precursor[3].ci)
    print("c4",chinRDynamic1D[this.moduleIndex].Precursor[4].ci)
    print("c5",chinRDynamic1D[this.moduleIndex].Precursor[5].ci)
    print("c6",chinRDynamic1D[this.moduleIndex].Precursor[6].ci)
end

function ReactorDynamics1DClass.Step(this,dt)
    chiReactorDynamics1DStep(this.moduleIndex,dt)
    chiReactorDynamics1DQuery(4,this.moduleIndex);
    --print("Power",chinRDynamic1D[this.moduleIndex].initPower)
    --print("Rho",chinRDynamic1D[this.moduleIndex].initRho)
    --print("c1",chinRDynamic1D[this.moduleIndex].Precursor[1].ci)
    --print("c2",chinRDynamic1D[this.moduleIndex].Precursor[2].ci)
    --print("c3",chinRDynamic1D[this.moduleIndex].Precursor[3].ci)
    --print("c4",chinRDynamic1D[this.moduleIndex].Precursor[4].ci)
    --print("c5",chinRDynamic1D[this.moduleIndex].Precursor[5].ci)
    --print("c6",chinRDynamic1D[this.moduleIndex].Precursor[6].ci)
end