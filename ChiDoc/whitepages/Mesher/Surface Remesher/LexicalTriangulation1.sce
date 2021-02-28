clear all;
clc;
funcprot(0)
rand('seed',200)

getd()

printf("Hello\n")

//========================================================== Generated points
N=10;               //Number of points
points=rand(N,2);
if (%F) then

points(1,1)  =0; points(1,2)  =0;
points(2,1)  =1; points(2,2)  =0;

points(3,1)  =0.25; points(3,2)  =0.25;
points(4,1)  =0.50; points(4,2)  =0.25;
points(5,1)  =0.50; points(5,2)  =0.50;
points(6,1)  =0.25; points(6,2)  =0.50;

points(N-1,1)=1; points(N-1,2)=1;
points(N,1)  =0; points(N,2)  =1;
end
matrixTime=500;
matrixPlot=%T;

//========================================================== Lexicographically sort points
[lexlist,lexpoints]=SortLexicographically2D(points);

//========================================================== Create initial convex hull
//The convex hull data structure is as follows:
// point_a point_b point_c neigbor_ab neighbor_bc neighbor_ca
//"point_" is only an index to a vertex not an actual xyz value
//"neighbor_" is the index of a triangle in the convexHull list
hullPoint1 =  lexlist(1);
hullPoint2 = -1;
hullPoint3 = -1;
for i=2:N
    if (points(lexlist(1),1)~=points(lexlist(i),1)) then
        hullPoint2=lexlist(i);
        break;
    end
end
for i=2:N
    a=points(lexlist(i),:);
    b=points(hullPoint1,:);
    c=points(hullPoint2,:);
    orientation = Orient2D(a,b,c);
    disp(orientation)
    
    if (   (lexlist(i)~=hullPoint2)   )
        hullPoint3=lexlist(i);
        if (orientation>0) then
            tempHullpoint=hullPoint2;
            hullPoint2 = hullPoint3;
            hullPoint3 = tempHullpoint;
            disp("Vertices flipped")
        end
        
        break;
    end
end
convexHull=[hullPoint1 hullPoint2 hullPoint3 -1 -1 -1]
[cH_itemcount cH_xsize] = size(convexHull)

scf(1)
clf(1)
if (matrixPlot)
     //****************************************************
     for t=1:(size(convexHull)(1))
         firstTri=[
         points(convexHull(t,1),:)
         points(convexHull(t,2),:)
         points(convexHull(t,3),:)
         points(convexHull(t,1),:)
         ]
         plot2d(firstTri(:,1),firstTri(:,2))
     end 
     scatter(lexpoints(1:N,1),lexpoints(1:N,2),,"black",".")
     a=gca();
     a.axes_visible = ["off" "off" "off"];
     a.box = "on"
     a.data_bounds = [-0.1,-0.1;1.1,1.1]
     
     sleep(matrixTime)
     //***************************************************
end


//========================================================== Create the left over point list
[unusedLexlist]=GetUnusedVertices(lexlist,convexHull)
hullEdges=FindOpenEdges(convexHull);

//========================================================== Iterate until all vertices are used
//while (unusedCount>0)

stopLoop=%F
iter=0;
while (~stopLoop)

    //=============== Convexify hull
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
            a=points(edgeCombo(2,2),:);        
            b=points(edgeCombo(1,1),:);
            c=points(edgeCombo(1,2),:);
            
            orientation = Orient2D(a,b,c);
            
            if (orientation>0) then
                concaveEdgesFound=%T;
                disp("Found concave edges")
                //Tell the owner of the edge that it now has a neighbor
                ownerTriangle1 = edgeCombo(1,3);
                ownerEdgeNo1   = edgeCombo(1,4);
                
                ownerTriangle2 = edgeCombo(2,3);
                ownerEdgeNo2   = edgeCombo(2,4);
                
                convexHull(ownerTriangle1,3+ownerEdgeNo1)=size(convexHull)(1)+1;
                convexHull(ownerTriangle2,3+ownerEdgeNo2)=size(convexHull)(1)+1;
                
                //Tell the new triangle which neighbor its attaching to
                newTriangle = [edgeCombo(1,1) edgeCombo(2,2) edgeCombo(1,2) -1 ownerTriangle2 ownerTriangle1]
                
                convexHull= [convexHull; newTriangle]; //Add triangle
                printf("Triangle added, number of unused verts=%d\n",size(unusedLexlist)(1))  
                
                //Update unused vertices and openEdges
                [unusedLexlist]=GetUnusedVertices(lexlist,convexHull)
                hullEdges=FindOpenEdges(convexHull);
                
               
               if (matrixPlot)
                    //****************************************************
                    for t=1:(size(convexHull)(1))
                        firstTri=[
                        points(convexHull(t,1),:)
                        points(convexHull(t,2),:)
                        points(convexHull(t,3),:)
                        points(convexHull(t,1),:)
                        ]
                        plot2d(firstTri(:,1),firstTri(:,2))
                    end 
                    scatter(lexpoints(1:N,1),lexpoints(1:N,2),,"black",".")
                    a=gca();
                    a.axes_visible = ["off" "off" "off"];
                    a.box = "on"
                    a.data_bounds = [-0.1,-0.1;1.1,1.1]
                    
                    sleep(matrixTime)
                    //***************************************************
                end
                break;
            end                  
        end
    end
    
    

    //=============== run over unused verts attempt to attach them to hull
    //for i=1:(size(unusedLexlist)(1))
    i=1;
    for t=1:(size(hullEdges)(1))
        a=points(unusedLexlist(i),:);
        b=points(hullEdges(t,1),:);
        c=points(hullEdges(t,2),:);
        
        orientation = Orient2D(a,b,c);
        
        if (orientation>0) then
            //Tell the owner of the edge that it now has a neighbor
            ownerTriangle = hullEdges(t,3);
            ownerEdgeNo   = hullEdges(t,4);
            convexHull(ownerTriangle,3+ownerEdgeNo)=size(convexHull)(1)+1;
            
            //Tell the new triangle which neighbor its attaching to
            newTriangle = [hullEdges(t,1) unusedLexlist(i) hullEdges(t,2) -1 -1 ownerTriangle]
            
            convexHull= [convexHull; newTriangle]; //Add triangle
            printf("Triangle added, number of unused verts=%d\n",size(unusedLexlist)(1))  
            
            //Update unused vertices and openEdges
            [unusedLexlist]=GetUnusedVertices(lexlist,convexHull)
            hullEdges=FindOpenEdges(convexHull);
           
            if (matrixPlot)
                //****************************************************
                for t=1:(size(convexHull)(1))
                    firstTri=[
                    points(convexHull(t,1),:)
                    points(convexHull(t,2),:)
                    points(convexHull(t,3),:)
                    points(convexHull(t,1),:)
                    ]
                    plot2d(firstTri(:,1),firstTri(:,2))
                end 
                scatter(lexpoints(1:N,1),lexpoints(1:N,2),,"black",".")
                a=gca();
                a.axes_visible = ["off" "off" "off"];
                a.box = "on"
                a.data_bounds = [-0.1,-0.1;1.1,1.1]
                //***************************************************
                sleep(matrixTime)
            end
            break;
        end                      
        if ( (size(unusedLexlist)(1)==0) )
            break;
        end
    end

    
    iter=iter+1;
    if ( (size(unusedLexlist)(1)==0) | (iter>2000))
        stopLoop=%T;
    end
    
    
end
disp("Iterations finished")

//=============== Convexify hull
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
        a=points(edgeCombo(2,2),:);        
        b=points(edgeCombo(1,1),:);
        c=points(edgeCombo(1,2),:);
        
        orientation = Orient2D(a,b,c);
        
        if (orientation>0) then
            concaveEdgesFound=%T;
            disp("Found concave edges")
            //Tell the owner of the edge that it now has a neighbor
            ownerTriangle1 = edgeCombo(1,3);
            ownerEdgeNo1   = edgeCombo(1,4);
            
            ownerTriangle2 = edgeCombo(2,3);
            ownerEdgeNo2   = edgeCombo(2,4);
            
            convexHull(ownerTriangle1,3+ownerEdgeNo1)=size(convexHull)(1)+1;
            convexHull(ownerTriangle2,3+ownerEdgeNo2)=size(convexHull)(1)+1;
            
            //Tell the new triangle which neighbor its attaching to
            newTriangle = [edgeCombo(1,1) edgeCombo(2,2) edgeCombo(1,2) -1 ownerTriangle2 ownerTriangle1]
            
            convexHull= [convexHull; newTriangle]; //Add triangle
            printf("Triangle added, number of unused verts=%d\n",size(unusedLexlist)(1))  
            
            //Update unused vertices and openEdges
            [unusedLexlist]=GetUnusedVertices(lexlist,convexHull)
            hullEdges=FindOpenEdges(convexHull);
            
            if (matrixPlot)
                //****************************************************
                for t=1:(size(convexHull)(1))
                    firstTri=[
                    points(convexHull(t,1),:)
                    points(convexHull(t,2),:)
                    points(convexHull(t,3),:)
                    points(convexHull(t,1),:)
                    ]
                    plot2d(firstTri(:,1),firstTri(:,2))
                end 
                scatter(lexpoints(1:N,1),lexpoints(1:N,2),,"black",".")
                a=gca();
                a.axes_visible = ["off" "off" "off"];
                a.box = "on"
                a.data_bounds = [-0.1,-0.1;1.1,1.1]
                //***************************************************
                sleep(matrixTime)
            end
            break;
        end                  
    end
end




scf(0)
clf(0)

subplot(221)
scatter(points(1:N,1),points(1:N,2),,"black",".")
a=gca();
a.axes_visible = ["off" "off" "off"];
//a.box = "off"
a.data_bounds = [-0.1,-0.1;1.1,1.1]

subplot(222)
plot2d(points(1:N,1),points(1:N,2))
scatter(points(1:N,1),points(1:N,2),,"black",".")
a=gca();
a.axes_visible = ["off" "off" "off"];
a.box = "on"
a.data_bounds = [-0.1,-0.1;1.1,1.1]


subplot(223)
plot2d(lexpoints(1:N,1),lexpoints(1:N,2))
scatter(lexpoints(1:N,1),lexpoints(1:N,2),,"black",".")
a=gca();
a.axes_visible = ["off" "off" "off"];
a.box = "on"
a.data_bounds = [-0.1,-0.1;1.1,1.1]


subplot(224)
for t=1:(size(convexHull)(1))
    firstTri=[
    points(convexHull(t,1),:)
    points(convexHull(t,2),:)
    points(convexHull(t,3),:)
    points(convexHull(t,1),:)
    ]
    plot2d(firstTri(:,1),firstTri(:,2))
end
    
    
scatter(lexpoints(1:N,1),lexpoints(1:N,2),,"black",".")
a=gca();
a.axes_visible = ["off" "off" "off"];
a.box = "on"
a.data_bounds = [-0.1,-0.1;1.1,1.1]

printf("Bye\n")



