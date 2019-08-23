
chiBindScene(0);
chiSetWindowProperties(chinGlobal.dwindowxsize,chinGlobal.dwindowysize,chinGlobal.dwindowxposi,chinGlobal.dwindowyposi);
chinGlobal.dwindowxsize,chinGlobal.dwindowysize=chiGetWindowProperties();
chiWindowCreate("CHI_TECHD");
chiBindScene(1);
chiSetWindowProperties("CHILD_OF",0);
chiBindScene(0);



