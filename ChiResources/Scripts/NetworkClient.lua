clientOn = false;
function ClientCommand()

serverIP = "127.0.0.1";

    if (WM_MBUTTONDOWN.occured) then
        print("Attempting to connect to server at "..serverIP.."...");
        result = chiNetworkStartClient(serverIP);
        if (result) then
            print("Success");
            clientOn = true;
        else
            print("Failure");
        end
    end
    
    if (clientOn) then
        --chiNetworkReceiveMessage();
        chiNetworkReceivePacket();
    end
end