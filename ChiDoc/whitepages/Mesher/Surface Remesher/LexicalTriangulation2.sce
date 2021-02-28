clear all;
clc;
funcprot(0)
//rand('seed',200)
rand('seed',200)

getd()

printf("Hello\n")

//========================================================== Generated points
N=10;               //Number of points
points=rand(N,2);

if (%F) then
    N=8;
    points=rand(N,2);
    points(1,1)  =0; points(1,2)  =0;
    points(2,1)  =1; points(2,2)  =0;
    del=0.25;
    points(3,1)  =del; points(3,2)  =del;
    points(4,1)  =1-del-del; points(4,2)  =del;
    points(5,1)  =1-del-del; points(5,2)  =1-del-del;
    points(6,1)  =del; points(6,2)  =1-del-del;
    
    points(N-1,1)=1; points(N-1,2)=1;
    points(N,1)  =0; points(N,2)  =1;
end
    

//========================================================== Lexicographically sort points
printf("================ Lexicographically sorting points\n");
[lexlist,lexpoints]=SortLexicographically2D(points);


//========================================================== Create initial convex hull
//The convex hull data structure is as follows:
// point_a point_b point_c neigbor_ab neighbor_bc neighbor_ca
//"point_" is only an index to a vertex not an actual xyz value
//"neighbor_" is the index of a triangle in the convexHull list
printf("================ Creating initial convex hull (i.e. first Triangle)\n");
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


//========================================================== Create the left over point list
[unusedLexlist]=GetUnusedVertices(lexlist,convexHull)

//========================================================== Iterate until all vertices are used
printf("================ Iterating over unused vertices\n");
stopLoop=%F
iter=0;
while (~stopLoop)
    //============================================ Convexify hull
    convexHull = ConvexifyHull(convexHull,points);
    
    //============================================ Run over unused verts attempt to attach them to hull
    [convexHull,unusedLexlist]=AttachUnusedVertices(convexHull,lexlist,unusedLexlist,points)

    iter=iter+1;
    if ( (size(unusedLexlist)(1)==0) | (iter>2000))
        stopLoop=%T;
    end
end
printf("================ Done iterating, final convexifying\n");

//========================================================== Final Convexify hull
convexHull = ConvexifyHull(convexHull,points);
lexi_convexHull=convexHull;
printf("================ Done convexifying the hull\n");
ExportAsOBJ("TestSurface4.obj", points,convexHull);

//========================================================== Create list of non-locally-delaunay edges
non_loc_del_edges = ListNonLocallyDelaunayEdges(convexHull,points);
temp_nlde=non_loc_del_edges //for plotting
temp_convexHull=convexHull

//========================================================== Iterate to remove non-locally delaunay edges
iter=0;
while (size(non_loc_del_edges)(1)>0)
//while (%F)
    iter=iter+1;
    printf("============ ITERATION %3d ==============\n",iter)
    [convexHull]=EdgeFlip(convexHull,non_loc_del_edges)
    non_loc_del_edges = ListNonLocallyDelaunayEdges(convexHull,points);
    if (iter==0) then
        //temp_nlde=non_loc_del_edges
        //temp_convexHull=convexHull
        disp(non_loc_del_edges)
        disp(convexHull)
    end
    
end





//temp_nlde=non_loc_del_edges
//while (size(non_loc_del_edges)(1)>0)
    
//end



scf(0)
clf(0)

subplot(321)
scatter(points(1:N,1),points(1:N,2),,"black",".")
a=gca();
a.axes_visible = ["off" "off" "off"];
//a.box = "off"
a.data_bounds = [-0.1,-0.1;1.1,1.1]

subplot(322)
plot2d(points(1:N,1),points(1:N,2))
scatter(points(1:N,1),points(1:N,2),,"black",".")
a=gca();
a.axes_visible = ["off" "off" "off"];
a.box = "on"
a.data_bounds = [-0.1,-0.1;1.1,1.1]


subplot(323)
plot2d(lexpoints(1:N,1),lexpoints(1:N,2))
scatter(lexpoints(1:N,1),lexpoints(1:N,2),,"black",".")
a=gca();
a.axes_visible = ["off" "off" "off"];
a.box = "on"
a.data_bounds = [-0.1,-0.1;1.1,1.1]


subplot(324)
for t=1:(size(lexi_convexHull)(1))
    firstTri=[
    points(lexi_convexHull(t,1),:)
    points(lexi_convexHull(t,2),:)
    points(lexi_convexHull(t,3),:)
    points(lexi_convexHull(t,1),:)
    ]
    plot2d(firstTri(:,1),firstTri(:,2))
end
    
    
scatter(lexpoints(1:N,1),lexpoints(1:N,2),,"black",".")
a=gca();
a.axes_visible = ["off" "off" "off"];
a.box = "on"
a.data_bounds = [-0.1,-0.1;1.1,1.1]


subplot(325)
for t=1:(size(temp_convexHull)(1))
    firstTri=[
    points(temp_convexHull(t,1),:)
    points(temp_convexHull(t,2),:)
    points(temp_convexHull(t,3),:)
    points(temp_convexHull(t,1),:)
    ]
    id=color("green")
    plot2d(firstTri(:,1),firstTri(:,2),id)
    xstring(mean(firstTri(:,1)),mean(firstTri(:,2)),string(t))
     id=color("gray")
    t=get("hdl")
    t.font_foreground=id
end

id=color("red")
for k=1:size(temp_nlde)(1)
    edge=temp_nlde(k,:)
    plotPoints=[points(edge(1),:); points(edge(2),:)];
    plot2d(plotPoints(:,1),plotPoints(:,2),id)
end

scatter(lexpoints(1:N,1),lexpoints(1:N,2),,"black",".")
dx=0.00
dy=-0.05
xstring(lexpoints(1:N,1)+dx,lexpoints(1:N,2)+dy,string(lexlist(1:N)))
a=gca();
a.axes_visible = ["off" "off" "off"];
a.box = "on"
a.data_bounds = [-0.1,-0.1;1.1,1.1]


subplot(326)
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






scf(1)
clf(1)
for t=1:(size(convexHull)(1))
    firstTri=[
    points(convexHull(t,1),:)
    points(convexHull(t,2),:)
    points(convexHull(t,3),:)
    points(convexHull(t,1),:)
    ]
    //id=color("green")
    plot2d(firstTri(:,1),firstTri(:,2))
    dx=-0.005
    dy=-0.015
    xstring(mean(firstTri(:,1))+dx,mean(firstTri(:,2))+dy,string(t))
    id=color("gray")
    t=get("hdl")
    t.font_foreground=id
end

circsize = 4;
xc=zeros(circsize,1); yc=zeros(circsize,1); r=zeros(circsize,1);
ex= 0.588; ey= 0.695;
xc(1)= 0.664; yc(1)= 0.662; r(1)= 0.163;
xc(2)= 0.571; yc(2)= 0.692; r(2)= 0.196;
xc(3)= 0.382; yc(3)= 0.555; r(3)= 0.282;
xc(4)= 0.382; yc(4)= 0.656; r(4)= 0.286;


scatter(lexpoints(1:N,1),lexpoints(1:N,2),,"black",".")

scatter(xc,yc,,"red","+")
scatter([ex; ex],[ey; ey],,"green",".")
dx=0.005
dy=-0.01
xstring(lexpoints(1:N,1)+dx,lexpoints(1:N,2)+dy,string(lexlist(1:N)-1))
a=gca();
//a.axes_visible = ["off" "off" "off"];
a.box = "on"
a.data_bounds = [-0.1,-0.1;1.1,1.1]

scf(2)
clf(2)

ponts = zeros(10,3);
ponts(1,1)= 0.018; ponts(1,2)= 0.977; ponts(1,3)=1
ponts(2,1)= 0.140; ponts(2,2)= 0.699; ponts(2,3)=2
ponts(3,1)= 0.515; ponts(3,2)= 0.306; ponts(3,3)=3
ponts(4,1)= 0.907; ponts(4,2)= 0.976; ponts(4,3)=4
ponts(5,1)= 0.551; ponts(5,2)= 0.887; ponts(5,3)=5
ponts(6,1)= 0.835; ponts(6,2)= 0.043; ponts(6,3)=6
ponts(7,1)= 0.751; ponts(7,2)= 0.279; ponts(7,3)=7
ponts(8,1)= 0.625; ponts(8,2)= 0.504; ponts(8,3)=8
ponts(9,1)= 0.726; ponts(9,2)= 0.812; ponts(9,3)=9
ponts(10,1)= 0.771; ponts(10,2)= 0.539; ponts(10,3)=10
cHull = [ 1 4 0];
cHull=[cHull; 1 7 4]
cHull=[cHull; 7 1 2]
cHull=[cHull; 7 8 4]
cHull=[cHull; 2 6 7]
cHull=[cHull; 7 9 8]
cHull=[cHull; 9 7 6]
cHull=[cHull; 2 5 6]
cHull=[cHull; 5 9 6]
cHull=[cHull; 5 3 9]
cHull=[cHull; 3 8 9]
cHull=[cHull; 3 4 8]
cHull=[cHull; 3 0 4]

for t=1:(size(cHull)(1))
    firstTri=[
    ponts(cHull(t,1)+1,:)
    ponts(cHull(t,2)+1,:)
    ponts(cHull(t,3)+1,:)
    ponts(cHull(t,1)+1,:)
    ]
    //id=color("green")
    plot2d(firstTri(:,1),firstTri(:,2))
    dx=-0.005
    dy=-0.015
    xstring(mean(firstTri(:,1))+dx,mean(firstTri(:,2))+dy,string(t))
    id=color("gray")
    t=get("hdl")
    t.font_foreground=id
end

scatter(ponts(1:N,1),ponts(1:N,2),,"black",".")
dx=0.005
dy=-0.01
xstring(ponts(1:N,1)+dx,ponts(1:N,2)+dy,string(ponts(1:N,3)-1))
a=gca();
//a.axes_visible = ["off" "off" "off"];
a.box = "on"
a.data_bounds = [-0.1,-0.1;1.1,1.1]


printf("Bye\n")



