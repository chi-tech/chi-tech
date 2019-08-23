serverOn = false;
function ServerCommand()

    if (not serverOn) then
        if (WM_LBUTTONDOWN.occured) then
            print("Starting Server...");
            result = chiNetworkStartServer();
            if (result) then
                serverOn = true;
                print("Success");
                chiNetworkBroadcast("Welcome to the Funhouse");
            else
                print("Failure");
            end
        end
    end
    
    if (WM_RBUTTONDOWN.occured) then
        print("Attempting to Update Server...");
        result  = chiNetworkUpdateServer();
    end
    
    if (serverOn) then
        --chiNetworkReceiveMessage();
        chiNetworkReceivePacket();
    end
end