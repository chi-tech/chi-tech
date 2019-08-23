directXDeviceCount=0;
prevdirectXDeviceCount=0;
directXDeviceCount=chiDirectXGetDeviceCount()
function UpdateDXDevices()
    if (not (directXDeviceCount==prevdirectXDeviceCount)) then
        prevdirectXDeviceCount=directXDeviceCount;
        
        chiDirectXGetList();
        for k=0,(directXDeviceCount-1) do
            output=string.format("Device %d, %s\n",k,chinDirectXDevices[k].name);
            output=output..string.format("Axes   : %d\n",chinDirectXDevices[k].axisCount);
            output=output..string.format("POVs   : %d\n",chinDirectXDevices[k].povCount);
            output=output..string.format("Buttons: %d  ",chinDirectXDevices[k].buttonCount);
            print(output)
        end
    end
    chiDirectXGetDeviceState(0);
    --print(chinDirectXDevices[0].button[0],chinDirectXDevices[0].button[1])
end