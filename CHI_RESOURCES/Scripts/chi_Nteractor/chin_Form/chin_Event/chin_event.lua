--Event pump
feature={}
feature.count=0;


--=============================================== Register feature
function dRegisterFeature(thefeature)
    feature.count=feature.count+1;
    k=feature.count;
    feature[k]=thefeature
end


--=============================================== Process Event Stack
cycleNumber=1;
function d_events()
    cycleNumber=cycleNumber+1;
    if (cycleNumber==1) then
        --chiSetWindowProperties("MINIMIZED");
        --chiSetWindowProperties("RESTORED");
    end
    
    if (WM_MOUSEMOVE.occured) then
        chiRequestSceneRefresh();
    end
    
    if (WM_SIZE.occured) then
        chinGlobal.dwindowxsize,chinGlobal.dwindowysize=chiGetWindowProperties();
        chiGraphicsCameraOrthoWidth(chinGlobal.dwindowxsize);
        chiGraphicsPositionCamera(chinGlobal.dwindowxsize/2,chinGlobal.dwindowysize/2,100.0);
    end

    if (WM_CHAR.occured) then
        if(WM_CHAR.iPar0 == 114) then
            --chiReloadShader("Picking"   );
            --chiReloadShader("Flat"      );
            --chiReloadShader("Phong"     );
            --chiReloadShader("Shadow1"   );
            --chiReloadShader("Shadow2"   );
            --chiReloadShader("PhongNS"   );
            --chiReloadShader("PhongNSNT" );
        end
    end

    for k=1,feature.count do
        if (feature[k]~=nil) then
            feature[k].ProcessEvents(feature[k]);
        end
    end
end

--=============================================== FilterEvents
function d_filter()
    --if ((WM_LBUTTONDOWN.occured  ) and (WM_LBUTTONDOWN.iPar4==0))   then WM_LBUTTONDOWN.occured  =false;  end
	--if ((WM_RBUTTONDOWN.occured  ) and (WM_RBUTTONDOWN.iPar4==0))   then WM_RBUTTONDOWN.occured  =false;  end
	--if ((WM_MBUTTONDOWN.occured  ) and (WM_MBUTTONDOWN.iPar4==0))   then WM_MBUTTONDOWN.occured  =false;  end
	--if ((WM_LBUTTONUP.occured    ) and (WM_LBUTTONUP.iPar4==0))     then WM_LBUTTONUP.occured    =false;  end
	--if ((WM_RBUTTONUP.occured    ) and (WM_RBUTTONUP.iPar4==0))     then WM_RBUTTONUP.occured    =false;  end
	--if ((WM_MBUTTONUP.occured    ) and (WM_MBUTTONUP.iPar4==0))     then WM_MBUTTONUP.occured    =false;  end
	--if ((WM_MOUSEWHEEL.occured   ) and (WM_MOUSEWHEEL.iPar4==0))    then WM_MOUSEWHEEL.occured   =false;  end
	--if ((WM_MOUSEMOVE.occured    ) and (WM_MOUSEMOVE.iPar4==0))     then WM_MOUSEMOVE.occured    =false;  end
	--if ((WM_SHIFTBUTTONDN.occured) and (WM_SHIFTBUTTONDN.iPar4==0)) then WM_SHIFTBUTTONDN.occured=false;  end
	--if ((WM_SHIFTBUTTONUP.occured) and (WM_SHIFTBUTTONUP.iPar4==0)) then WM_SHIFTBUTTONUP.occured=false;  end
	--if ((WM_CTLBUTTONDN.occured  ) and (WM_CTLBUTTONDN.iPar4==0))   then WM_CTLBUTTONDN.occured  =false;  end
	--if ((WM_CTLBUTTONUP.occured  ) and (WM_CTLBUTTONUP.iPar4==0))   then WM_CTLBUTTONUP.occured  =false;  end
	----if ((WM_KEYDN.occured        ) and (WM_KEYDN.iPar4==0))         then WM_KEYDN.occured        =false;  end
	----if ((WM_KEYUP.occured        ) and (WM_KEYUP.iPar4==0))         then WM_KEYUP.occured        =false;  end
	--if ((WM_MOUSELEAVE.occured   ) and (WM_MOUSELEAVE.iPar4==0))    then WM_MOUSELEAVE.occured   =false;  end
	--if ((WM_ENTERDN.occured      ) and (WM_ENTERDN.iPar4==0))       then WM_ENTERDN.occured      =false;  end
	--if ((WM_ENTERUP.occured      ) and (WM_ENTERUP.iPar4==0))       then WM_ENTERUP.occured      =false;  end
    --if ((WM_SIZE.occured         ) and (WM_ENTERUP.iPar4==0))       then WM_SIZE.occured         =false;  end
    
end