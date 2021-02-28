dofile("CHI_RESOURCES/Scripts/chi_Nteractor/chin_GlobalFunctions/chin_globalfs_getscriptpath.lua")
chiplot_baseDir = chinGetScriptPath();

chiplot_MAXPLOTS = 20;

chiplot = {}
for k=1,chiplot_MAXPLOTS do
    chiplot[k] = {};
    chiplot[k].initialized      = false;
    chiplot[k].dimension        = 2;
    chiplot[k].xlabel           = "x";
    chiplot[k].ylabel           = "y";
    chiplot[k].zlabel           = "z";
    
    chiplot[k].xautorange       = true;
    chiplot[k].yautorange       = true;
    chiplot[k].zautorange       = true;
    
    chiplot[k].xmajordiv        = 10;
    chiplot[k].ymajordiv        = 10;
    chiplot[k].zmajordiv        = 10;
    
    chiplot[k].xminordiv        = 2;
    chiplot[k].yminordiv        = 2;
    chiplot[k].zminordiv        = 2;
    
    chiplot[k].xmax             = 1.0;
    chiplot[k].ymax             = 1.0;
    chiplot[k].zmax             = 1.0;
    
    chiplot[k].xmin             = 0.0;
    chiplot[k].ymin             = 0.0;
    chiplot[k].zmin             = 0.0;
    
    chiplot[k].xislog           = false;
    chiplot[k].yislog           = false;
    chiplot[k].zislog           = false;
end


--######################################################### Plot array
function chiplotArray(theArray,pnum)
    if (chiplot[k].initialized == false) then

    end

end