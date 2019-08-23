--What are panels? Panels are essentially viewports. Stuff drawn within panels 
--are clipped on the edges. Panels can be free floating or bound by edges.
--Panels should have different edge styles and should always have a resize 
--handler.
--Panels can also be constrained to a minimum size.
--For connected panels, code should be implemented that dynamically moves all
--panels that are attached together.

PanelClass = {}
PanelClass.__index = PanelClass
PanelCount = 0;

PanelSurface = chiLoadSurface(chinPanelDir .. "chin_panel.obj");
PanelTexture = chiLoadTexture("CHI_RESOURCES/Textures/chin_Icons.png");


dofile(chinPanelDir .. "chin_panel_00_constrdestr.lua")
dofile(chinPanelDir .. "chin_panel_01_functions.lua")