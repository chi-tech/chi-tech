chiBindScene(0);
chiSetWindowProperties(chinGlobal.dwindowxsize,chinGlobal.dwindowysize,dwindowxposi,dwindowyposi);
chinGlobal.dwindowxsize,chinGlobal.dwindowysize = chiGetWindowProperties();
--chiWindowCreate("chi_Nplot");
--chiBindScene(1);
--chiSetWindowProperties("CHILD_OF",0);
chiBindScene(0);

chiLightCreate("DefaultLight");
ambient=0.0;
chiLightSetProperty("DefaultLight","Ambient",ambient,ambient,ambient)
chiLightSetProperty("DefaultLight","Diffuse",ambient,ambient,ambient)
chiLightSetProperty("DefaultLight","Attenuation",1000.0,0.001,0.0)
chiLightSetProperty("DefaultLight","Position",0.0,-280,138)


