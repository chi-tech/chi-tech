clear all;
clc;
funcprot(0)
printf("Hello\n")

//##############################################################################
function [rlexlist, tempPointList]=SortLexicographically2D(pointList)
    listSize = size(pointList)(1);
    tempPointList = pointList;
    rlexlist = linspace(1,listSize,listSize)
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

//##############################################################################
function valid = CheckValidTri2D(v1,v2,v3)
    vert1=[v1(1) v1(2) 0];
    vert2=[v2(1) v2(2) 0];
    vert3=[v3(1) v3(2) 0];
    
    AB = vert2-vert1;
    nAB = AB/norm(AB);
    
    BC = vert3-vert2;
    nBC = BC/norm(BC);
    
    if (abs(nAB*nBC')>0.999999) then
        valid=%F;
        return;
    end
    crossP = cross(AB,BC)
    crossP = crossP/norm(crossP)
    disp(crossP)
    if (crossP(3)<0) then
        //valid=%F;
        //return;
    end
    valid=%T
endfunction


N = 10;             //Number of points

points=[
0.0982278226	0.333627197
0.7653228931	0.2347061365
0.2558645175	0.8037821481
0.3081272112	0.8671842781
0.8975283975	0.1435997586
0.2429157385	0.1533325249
0.2663349204	0.4709250691
0.3667177816	0.2964881263
0.3963756175	0.5602469026
0.4162725839	0.6990938702
]

[lexlist,sorted]=SortLexicographically2D(points)
disp(lexlist)
disp(sorted)

scf(0);
clf();
scatter(points(1:N),points(1:N,2),,"black",".")
a=gca()
a.axes_visible = ["off" "off" "off"];
a.box = "off"

sleep(1000);

for i=1:(N-2)
    j=i+1; k=1;
    iter=0;
    foundTriangle=%F
    while (~foundTriangle)
        combo = [lexlist(i) lexlist(j) lexlist(j+k)];
        disp(combo)
        foundTriangle = CheckValidTri2D(points(combo(1),1:2),points(combo(2),1:2),points(combo(3),1:2))
        if (~foundTriangle) then
            k=k+1;
        end
        if (k>N) then
            break;
        end
        
        iter=iter+1;
        if (iter>100)
            foundTriangle=%T;
        end
    
    end
    if foundTriangle then
        combo = [lexlist(i) lexlist(j) lexlist(j+k)];
        pointsToPlot=[
        points(combo(1),1:2)
        points(combo(2),1:2)
        points(combo(3),1:2)
        points(combo(1),1:2)
        ]
        disp(pointsToPlot)
        plot2d(pointsToPlot(1:4,1),pointsToPlot(1:4,2))
    end
end
a=gca()
a.axes_visible = ["off" "off" "off"];
a.box = "off"


printf("Bye\n")



