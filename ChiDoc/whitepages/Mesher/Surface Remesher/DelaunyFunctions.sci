//##############################################################################
//Exports the convex hull as wavefront object
function []=ExportAsOBJ(fileName,in_points,in_convexHull)
    outputFile = mopen(fileName,"wt");
    mfprintf(outputFile,"# Exported mesh file from tringulation script\n");
    mfprintf(outputFile,"o ConvexHull\n");
    for p=1:size(in_points)(1)
        mfprintf(outputFile,"v %9.6f %9.6f  0.000000\n",in_points(p,1),in_points(p,2));
    end
    mfprintf(outputFile,"vn 0.0000 0.0000 1.0000\n");
    mfprintf(outputFile,"s off\n");
    
    for t=1:size(in_convexHull)(1)
        mfprintf(outputFile,"f %d//1 %d//1 %d//1\n",in_convexHull(t,1),in_convexHull(t,2),in_convexHull(t,3) );
    end
    
    mclose(outputFile);
    printf("FILE %s EXPORTED\n",fileName);
endfunction



//##############################################################################
//Lexicographically sort the vertices first by x then by y
//Returns rlexlist      A vector of the indices of the sorted list to the original list
//        tempPointList A sorted list of points (x,y)
function [rlexlist, tempPointList]=SortLexicographically2D(pointList)
    listSize = size(pointList)(1);
    tempPointList = pointList;
    rlexlist = linspace(1,listSize,listSize)';
    //========== Sort by x
    for i=listSize:-1:1
        for j=(i-1):-1:1
            if (tempPointList(i,1)<tempPointList(j,1))
                tempPointx = tempPointList(j,1);
                tempPointy = tempPointList(j,2);
                tempIndex = rlexlist(j);
                
                tempPointList(j,1)=tempPointList(i,1);
                tempPointList(j,2)=tempPointList(i,2);
                rlexlist(j) = rlexlist(i)
                
                tempPointList(i,1)=tempPointx;
                tempPointList(i,2)=tempPointy;
                rlexlist(i) = tempIndex
            end
        end
    end
    
    //========== Sort by y
    for i=listSize:-1:1
        for j=(i-1):-1:1
            if ((abs(tempPointList(i,1)-tempPointList(j,1) )<0.0000001) & (tempPointList(i,2)<tempPointList(j,2)))
                tempPointx = tempPointList(j,1);
                tempPointy = tempPointList(j,2);
                tempIndex = rlexlist(j);
                
                tempPointList(j,1)=tempPointList(i,1);
                tempPointList(j,2)=tempPointList(i,2);
                rlexlist(j) = rlexlist(i)
                
                tempPointList(i,1)=tempPointx;
                tempPointList(i,2)=tempPointy;
                rlexlist(i) = tempIndex
            end
        end
    end
endfunction

//############################################################################## Orient 2D
function result=Orient2D(a,b,c)
    ba=[(a-b) 0];
    bc=[(c-b) 0];
    ac_cross_bc = cross(ba,bc);
    result = ac_cross_bc(3);
endfunction

//############################################################################## Create the left over point list
function [out_unusedLexlist]=GetUnusedVertices(in_lexlist,in_triangles)
    [pcount,pcountx]=size(in_lexlist);
    [tcount,tcountx]=size(in_triangles);

    out_unusedLexlist   = [];
    for i=1:pcount
        foundInHull=%F
        for j=1:tcount
            if (   (in_lexlist(i)==in_triangles(j,1)) | (in_lexlist(i)==in_triangles(j,2)) | (in_lexlist(i)==in_triangles(j,3)) )
                foundInHull=%T;
            end
        end
        if (~foundInHull) then
            out_unusedLexlist=[out_unusedLexlist; in_lexlist(i)];
        end
    end
endfunction



//############################################################################## Find Open Edges
function [out_convexEdges]=FindOpenEdges(in_convexHull)
    temp_openEdges=[];
    for t=1:size(in_convexHull)(1)
        ownerTriangle=t;
        edges=[
                [in_convexHull(t,1) in_convexHull(t,2) in_convexHull(t,4)]
                [in_convexHull(t,2) in_convexHull(t,3) in_convexHull(t,5)]
                [in_convexHull(t,3) in_convexHull(t,1) in_convexHull(t,6)]   ];
        for e=1:3
            if (edges(e,3)<0) then //valid edge
                temp_openEdges=[temp_openEdges;edges(e,1:2) ownerTriangle e];
            end
        end
    end
    
    out_convexEdges = temp_openEdges(1,1:4);
    last_vertex=out_convexEdges(1,2);
    
    while (size(out_convexEdges)(1)~=  size(temp_openEdges)(1) )
        for e=1:size(temp_openEdges)(1)
            if (temp_openEdges(e,1)==last_vertex) then
                out_convexEdges=[out_convexEdges; temp_openEdges(e,1:4)];
                last_vertex = temp_openEdges(e,2);
            end
        end
    end
endfunction



//############################################################################## ConvexifyHull
function [out_convexHull]=ConvexifyHull(in_convexHull,in_points)
    out_convexHull = in_convexHull;
    hullEdges=FindOpenEdges(out_convexHull);
    
    concaveEdgesFound=%T;
    while (concaveEdgesFound)
        concaveEdgesFound=%F;
        for t=1:(size(hullEdges)(1))
            edgeCombo = [hullEdges(t,:)];
                         
            if (t<(size(hullEdges)(1))) then
                edgeCombo = [edgeCombo; hullEdges(t+1,:)];
            else
                edgeCombo = [edgeCombo; hullEdges(1,:)];
            end
            a=in_points(edgeCombo(2,2),:);        
            b=in_points(edgeCombo(1,1),:);
            c=in_points(edgeCombo(1,2),:);
            
            orientation = Orient2D(a,b,c);
            
            if (orientation>0) then
                concaveEdgesFound=%T;
                printf("Found concave edges\n")
                //Tell the owner of the edge that it now has a neighbor
                ownerTriangle1 = edgeCombo(1,3);
                ownerEdgeNo1   = edgeCombo(1,4);
                
                ownerTriangle2 = edgeCombo(2,3);
                ownerEdgeNo2   = edgeCombo(2,4);
                
                out_convexHull(ownerTriangle1,3+ownerEdgeNo1)=size(out_convexHull)(1)+1;
                out_convexHull(ownerTriangle2,3+ownerEdgeNo2)=size(out_convexHull)(1)+1;
                
                //Tell the new triangle which neighbor its attaching to
                newTriangle = [edgeCombo(1,1) edgeCombo(2,2) edgeCombo(1,2) -1 ownerTriangle2 ownerTriangle1]
                
                out_convexHull= [out_convexHull; newTriangle]; //Add triangle
                printf("Triangle added in ConvexifyHull\n")  
                
                //OpenEdges
                hullEdges=FindOpenEdges(out_convexHull);
                break;
            end                  
        end
    end
endfunction

//############################################################################## Attached unused vertices
function [out_convexHull,out_unusedLexlist]=AttachUnusedVertices(in_convexHull,in_lexlist,in_unusedLexlist,in_points)
    out_convexHull=in_convexHull;
    out_unusedLexlist=in_unusedLexlist;
    hullEdges=FindOpenEdges(out_convexHull);
    i=1;
    for t=1:(size(hullEdges)(1))
        a=points(out_unusedLexlist(i),:);
        b=points(hullEdges(t,1),:);
        c=points(hullEdges(t,2),:);
        
        orientation = Orient2D(a,b,c);
        
        if (orientation>0) then
            //Tell the owner of the edge that it now has a neighbor
            ownerTriangle = hullEdges(t,3);
            ownerEdgeNo   = hullEdges(t,4);
            out_convexHull(ownerTriangle,3+ownerEdgeNo)=size(out_convexHull)(1)+1;
            
            //Tell the new triangle which neighbor its attaching to
            newTriangle = [hullEdges(t,1) out_unusedLexlist(i) hullEdges(t,2) -1 -1 ownerTriangle]
            
            out_convexHull= [out_convexHull; newTriangle]; //Add triangle
            printf("Triangle added, number of unused verts=%d\n",size(out_unusedLexlist)(1)-1)  
            //Update unused vertices and openEdges
            [out_unusedLexlist]=GetUnusedVertices(in_lexlist,out_convexHull)
            
            hullEdges=FindOpenEdges(out_convexHull);
         
            break;
        end                      
        if ( (size(out_unusedLexlist)(1)==0) )
            break;
        end
    end
endfunction


//############################################################################## Perform in-circle test
//Returns >0 if d is inside the circum-circle of abc
//Returns =0 if d is on the circum-circle of abc
//Returns <0 if d is outside the circum-circle of abc
function [out_value]=InCircle(a,b,c,d)
    m = [a(1) a(2) a(1)^2+a(2)^2 1
         b(1) b(2) b(1)^2+b(2)^2 1
         c(1) c(2) c(1)^2+c(2)^2 1
         d(1) d(2) d(1)^2+d(2)^2 1];
    out_value = det(m);
endfunction



//############################################################################## Get list of non-locally delaunay edges
function [out_non_loc_del_edges] = ListNonLocallyDelaunayEdges(in_convexHull,in_points)
    out_non_loc_del_edges=[];    
    for t=1:(size(in_convexHull)(1))
        tau_1 = in_convexHull(t,:);
        for e=1:3
            if (tau_1(3+e)>-1) then   
                if (e<3) then
                    edge=[tau_1(e) tau_1(e+1) t tau_1(3+e)];
                else
                    edge=[tau_1(e) tau_1(1) t tau_1(3+e)];
                end
                
                tau_2_index = tau_1(3+e);
                if (tau_2_index>=0) then
                    tau_2 = in_convexHull(tau_2_index,:);
                    
                    a=in_points(tau_1(1),:);
                    b=in_points(tau_1(2),:);
                    c=in_points(tau_1(3),:);
                    
                    for v=1:3
                        if ( (tau_2(v)~=edge(1)) & (tau_2(v)~=edge(2))) then
                            d=in_points(tau_2(v),:);
                            break;
                        end
                    end
                    
                    localDelaunay = InCircle(a,b,c,d);
                    if (localDelaunay>0.000001) then
                        //============== Check for duplicates
                        duplicateFound=%F;
                        for d=1:size(out_non_loc_del_edges)(1)
                            curEdge=out_non_loc_del_edges(d,:);
                            if (  (curEdge(1)==edge(2))  &  (curEdge(2)==edge(1)) ) then
                                duplicateFound = %T;
                                break;
                            end
                        end
                        if (~duplicateFound) then
                            out_non_loc_del_edges=[out_non_loc_del_edges; edge];
                        end

                    end
                end
            end
        end
        
    end
    
endfunction


//############################################################################## Edge flip
function [out_convexHull]=EdgeFlip(in_convexHull,in_non_loc_del_edges)
    out_convexHull=in_convexHull;
    //for e=1:size(in_non_loc_del_edges)(1)
    for e=1:1
        //================================================== Get the edge's triangles
        curEdge=in_non_loc_del_edges(e,:);
        tau_1_index = curEdge(3);
        tau_2_index = curEdge(4);
        tau_1  =in_convexHull(curEdge(3),:);
        tau_2  =in_convexHull(curEdge(4),:);
        
        //================================================== Find which edge no this edge is for Tau 1 and 2
        tau_1_enum = 0;
        tau_2_enum = 0;
            if (   (  (curEdge(1)==tau_1(1)) & (curEdge(2)==tau_1(2))  ) | (  (curEdge(2)==tau_1(1)) & (curEdge(1)==tau_1(2))  )   ) then
                tau_1_enum = 1;
        elseif (   (  (curEdge(1)==tau_1(2)) & (curEdge(2)==tau_1(3))  ) | (  (curEdge(2)==tau_1(2)) & (curEdge(1)==tau_1(3))  )   ) then
                tau_1_enum = 2;
        elseif (   (  (curEdge(1)==tau_1(3)) & (curEdge(2)==tau_1(1))  ) | (  (curEdge(2)==tau_1(3)) & (curEdge(1)==tau_1(1))  )   ) then
                tau_1_enum = 3;
        end
            if (   (  (curEdge(1)==tau_2(1)) & (curEdge(2)==tau_2(2))  ) | (  (curEdge(2)==tau_2(1)) & (curEdge(1)==tau_2(2))  )   ) then
                tau_2_enum = 1;
        elseif (   (  (curEdge(1)==tau_2(2)) & (curEdge(2)==tau_2(3))  ) | (  (curEdge(2)==tau_2(2)) & (curEdge(1)==tau_2(3))  )   ) then
                tau_2_enum = 2;
        elseif (   (  (curEdge(1)==tau_2(3)) & (curEdge(2)==tau_2(1))  ) | (  (curEdge(2)==tau_2(3)) & (curEdge(1)==tau_2(1))  )   ) then
                tau_2_enum = 3;
        end
        
        //================================================== Find neighbour A and B
        neighbourA=-1;
        neighbourB=-1;
            if (tau_1_enum==1) then
            neighbourA = tau_1(3+2);
            neighbourB = tau_1(3+3);
        elseif (tau_1_enum==2) then
            neighbourA = tau_1(3+3);
            neighbourB = tau_1(3+1);
        elseif (tau_1_enum==3) then
            neighbourA = tau_1(3+1);
            neighbourB = tau_1(3+2);
        end
        
        //================================================== Find neighbour C and D
        neighbourC=-1;
        neighbourD=-1;
            if (tau_2_enum==1) then
            neighbourC = tau_2(3+2);
            neighbourD = tau_2(3+3);
        elseif (tau_2_enum==2) then
            neighbourC = tau_2(3+3);
            neighbourD = tau_2(3+1);
        elseif (tau_2_enum==3) then
            neighbourC = tau_2(3+1);
            neighbourD = tau_2(3+2);
        end

        //================================================== Find vertex A and B
        vertA=[];
        for v=1:3
            if (   (tau_1(v)~=curEdge(1))  &  (tau_1(v)~=curEdge(2))   ) then
                vertA=tau_1(v);
                break;
            end
        end
        vertB=[];
        for v=1:3
            if (   (tau_2(v)~=curEdge(1))  &  (tau_2(v)~=curEdge(2))   ) then
                vertB=tau_2(v);
                break;
            end
        end
        
        //================================================== Define modified triangles
        tau_A = [vertA vertB curEdge(2) tau_2_index neighbourD neighbourA]
        tau_B = [vertB vertA curEdge(1) tau_1_index neighbourB neighbourC]

        out_convexHull(curEdge(3),:)=tau_A(:);
        out_convexHull(curEdge(4),:)=tau_B(:);
        
        //================================================== Update neighbours
        for e=1:3
            if (neighbourB>0) then
                if (out_convexHull(neighbourB,3+e)==tau_1_index) then
                    out_convexHull(neighbourB,3+e)=tau_2_index;
                end
            end
            
            if (neighbourD>0) then
                if (out_convexHull(neighbourD,3+e)==tau_2_index) then
                    out_convexHull(neighbourD,3+e)=tau_1_index;
                end
            end
        end
        
        
    end
endfunction



