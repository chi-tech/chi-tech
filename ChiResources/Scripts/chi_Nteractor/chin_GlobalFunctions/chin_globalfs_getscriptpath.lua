--================================== Get current folder
--Obtains the folder of the currently exeucting
--script file.
function chinGetScriptPath()
   local filePath=debug.getinfo(2, "S").source:sub(2);
   local length=string.len(filePath);
   rfilePath=string.reverse(filePath);
   start=string.find(rfilePath,"\\");
   if (start==nil) then
        start=string.find(rfilePath,"/");
   end
   if (start==nil) then
        start=string.find(rfilePath,"//");
   end
   if (start==nil) then
        filePath=""
   end
   
   return string.sub(filePath,1,length-start);
end