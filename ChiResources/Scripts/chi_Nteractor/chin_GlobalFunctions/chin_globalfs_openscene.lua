function chinOpenScene(fileName)
    print("Opening Scene")
    chunk,err=assert(loadfile(fileName));
    print(err)
    if (chunk==nil) then print(err);
                else chunk() end
end