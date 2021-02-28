function SaveTables()
    io.write("\n--=========================== Tables\n")
    numTables = chiTableQuery(0);
    for tcount=1,numTables do
        k=tcount-1;
        chiTableQuery(4,k);
        io.write("temp = TableClass.New(\""..chinTable[k].name.."\");\n");
        io.write("chiTableLoadFromFile("..string.format("%d",k)..",\""..chinTable[k].fileName.."\");\n")
    end
    io.write("print(\"Tables loaded\");\n")
    
    
end