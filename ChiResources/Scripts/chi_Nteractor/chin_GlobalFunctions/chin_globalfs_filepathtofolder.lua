--================================== Get folder path from file path
--Obtains the folder of the indicated file.
function chinFilePathToFolder(filePath)
   local length=string.len(filePath);
   rfilePath=string.reverse(filePath);
   start=string.find(rfilePath,"\\");
   if (start==nil) then
        start=string.find(rfilePath,"/");
   end
   if (start==nil) then
        filePath=""
   end
   
   return string.sub(filePath,1,length-start);
end