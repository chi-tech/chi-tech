--######################################################### Process Events
function ManipulatorClass.ProcessEvents(this)
    currentScene=chiGetScene();
    chiBindScene(this.sceneScope);
    --======================= L mouse down
    if (WM_LBUTTONDOWN.occured) then
        this.lbuttondn = true;
        if     (WM_LBUTTONDOWN.iPar5==this.xobjNum)  then
            this.lbuttonscope = 1;
        elseif (WM_LBUTTONDOWN.iPar5==this.yobjNum)  then
            this.lbuttonscope = 2;
        elseif (WM_LBUTTONDOWN.iPar5==this.zobjNum)  then
            this.lbuttonscope = 3;
        else
            this.lbuttonscope = 0;
        end  
    end
    
    --======================= L mouse up
    if (WM_LBUTTONUP.occured) then
        this.lbuttondn = false;  
    end
    
    --======================= Mouse move
    if (WM_MOUSEMOVE.occured) then
        
        if (this.lbuttondn) then
            if (this.lbuttonscope==1) then
                x,y,z = chiTransformGet(this.manipTransform);
                
                mouse_dx = WM_MOUSEMOVE.iPar0 - WM_MOUSEMOVE.iPar2;
                mouse_dy = WM_MOUSEMOVE.iPar1 - WM_MOUSEMOVE.iPar3;
                
                --Window coords for unit movement

                wx1,wy1,wz1 = chi3DPointToScreen(x,y,z);
               
                wx2,wy2,wz2 = chi3DPointToScreen(x+1.0,y,z);
               
                
                sum1=       (wx2-wx1)*(wx2-wx1);
                sum1=sum1 + (wy2-wy1)*(wy2-wy1);
                sum1=sum1 + (wz2-wz1)*(wz2-wz1);
                Dr = math.sqrt(sum1); --[px/unit]
                
                Dw=math.sqrt(mouse_dx*mouse_dx + mouse_dy*mouse_dy);
                Dw=mouse_dx
                Du=Dw*2/Dr;
                
                if (wx2<wx1) then Du=Du*-1; end
               
                --chiTransformSetTranslation(this.manipTransform,x+Du, y,z);
                this.MasterTranslate(this,x+Du, y,z);
            end
            if (this.lbuttonscope==2) then
                x,y,z = chiTransformGet(this.manipTransform);
                
                mouse_dx = WM_MOUSEMOVE.iPar0 - WM_MOUSEMOVE.iPar2;
                mouse_dy = WM_MOUSEMOVE.iPar1 - WM_MOUSEMOVE.iPar3;
                
                --Window coords for unit movement

                wx1,wy1,wz1 = chi3DPointToScreen(x,y,z);
              
                wx2,wy2,wz2 = chi3DPointToScreen(x,y+1.0,z);
                
                
                sum1=       (wx2-wx1)*(wx2-wx1);
                sum1=sum1 + (wy2-wy1)*(wy2-wy1);
                sum1=sum1 + (wz2-wz1)*(wz2-wz1);
                Dr = math.sqrt(sum1); --[px/unit]
                
                Dw=math.sqrt(mouse_dx*mouse_dx + mouse_dy*mouse_dy);
                Dw=mouse_dx
                Du=Dw*2/Dr;
               if (wx2<wx1) then Du=Du*-1; end
                --chiTransformSetTranslation(this.manipTransform,x, y+Du,z);
                this.MasterTranslate(this,x, y+Du,z);
            end
            if (this.lbuttonscope==3) then
                x,y,z = chiTransformGet(this.manipTransform);
                
                mouse_dx = WM_MOUSEMOVE.iPar0 - WM_MOUSEMOVE.iPar2;
                mouse_dy = WM_MOUSEMOVE.iPar1 - WM_MOUSEMOVE.iPar3;
                
                --Window coords for unit movement

                wx1,wy1,wz1 = chi3DPointToScreen(x,y,z);

                wx2,wy2,wz2 = chi3DPointToScreen(x,y,z+1.0);

                
                sum1=       (wx2-wx1)*(wx2-wx1);
                sum1=sum1 + (wy2-wy1)*(wy2-wy1);
                sum1=sum1 + (wz2-wz1)*(wz2-wz1);
                Dr = math.sqrt(sum1); --[px/unit]
                
                Dw=math.sqrt(mouse_dx*mouse_dx + mouse_dy*mouse_dy);
                Dw=-mouse_dy
                Du=Dw*2/Dr;
               if (wy2<wy1) then Du=Du*-1; end
                --chiTransformSetTranslation(this.manipTransform,x, y,z+Du);
                this.MasterTranslate(this,x, y,z+Du);
            end
        end
    end
    
    --======================= Slave proc
    if (this.objectSlaveCount>0) then
        if (this.hidden) then
            this.UnHide(this);
        end
    end
    
    chiBindScene(currentScene);
end