-- |---------------------------------------------------------------|
-- | chi_Nplot Main Script                                         |
-- |---------------------------------------------------------------|

-- ========================================= Defining Directories
nBaseDir = "CHI_RESOURCES/Scripts/chi_Nplot/";
nCameDir = "CHI_RESOURCES/Scripts/chi_Nplot/chin_Camera/";
nWindDir = "CHI_RESOURCES/Scripts/chi_Nplot/chin_WindowSetup/";
nPaneDir = "CHI_RESOURCES/Scripts/chi_Nplot/chin_Panel/";
nEvenDir = "CHI_RESOURCES/Scripts/chi_Nplot/chin_Event/";



-- ========================================= Defining Includes
dofile(nBaseDir.."chi_Nplot.ini");
dofile(nWindDir.."chin_windowsetup.lua");
dofile(nEvenDir.."chin_event.lua");
dofile(nCameDir.."chin_camera.lua");
--dofile(nPaneDir.."chin_panel.lua");



-- ========================================= Defining Variables
xLabel,yLabel = chiGetWindowProperties();



-- ========================================= Defining Scenes

-- Creating Labels
chiSetLabel("x","I am the x-axis",xLabelx,xLabely,0,0,0,1,5);
chiSetLabel("y","I am the y-axis",yLabelx,yLabely,0,0,0,1,5);

-- Creating Axis
xAxis = chiLineCreate("X-Axis");
chiLineAddVertex(xAxis,xVertx,xVerty,0.0);
chiLineAddVertex(xAxis,xVertX,xVertY,0.0);

yAxis = chiLineCreate("Y-Axis");
chiLineAddVertex(yAxis,yVertx,yVerty,0.0);
chiLineAddVertex(yAxis,yVertX,yVertY,0.0);

-- Changing Line Properties
chiLineChangeColor(xAxis,0.0,0.0,0.0,1.0);
chiLineChangeColor(yAxis,0.0,0.0,0.0,1.0);

chiLineSetStipple(xAxis,false,0.0,2.0);
chiLineSetStipple(yAxis,false,0.0,2.0);


function chi_Nplot()
    chiBindScene(0);
    chinEvents();
    chinFilter();
end