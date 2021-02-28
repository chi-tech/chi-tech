clear all;
clc;
funcprot(0)
//rand('seed',200)
rand('seed',200)

scf(2)
clf(2)
ponts = zeros(8,3);
ponts(1,1)=-4.000; ponts(1,2)=-4.000; ponts(1,3)=1; ponts(1,4)= 1.000;
ponts(2,1)=-4.000; ponts(2,2)= 4.000; ponts(2,3)=2; ponts(2,4)= 1.000;
ponts(3,1)= 4.000; ponts(3,2)= 4.000; ponts(3,3)=3; ponts(3,4)= 1.000;
ponts(4,1)= 4.000; ponts(4,2)=-4.000; ponts(4,3)=4; ponts(4,4)= 1.000;
ponts(5,1)=-3.000; ponts(5,2)= 3.000; ponts(5,3)=5; ponts(5,4)= 1.000;
ponts(6,1)=-3.000; ponts(6,2)=-3.000; ponts(6,3)=6; ponts(6,4)= 1.000;
ponts(7,1)= 3.000; ponts(7,2)=-3.000; ponts(7,3)=7; ponts(7,4)= 1.000;
ponts(8,1)= 3.000; ponts(8,2)= 3.000; ponts(8,3)=8; ponts(8,4)= 1.000;
cHull = [ ]
cHull=[cHull; 3 6 0 2 1 -1 ]
cHull=[cHull; 6 5 0 -1 5 0 ]
cHull=[cHull; 3 7 6 6 -1 0 ]
cHull=[cHull; 1 1 1 1 1 1 ]
cHull=[cHull; 1 1 1 1 1 1 ]
cHull=[cHull; 4 0 5 9 1 -1 ]
cHull=[cHull; 3 2 7 -1 7 2 ]
cHull=[cHull; 2 4 7 8 -1 6 ]
cHull=[cHull; 2 1 4 -1 9 7 ]
cHull=[cHull; 1 0 4 -1 5 8 ]


loud = %T
for t=1:(size(cHull)(1))
    firstTri=[
    ponts(cHull(t,1)+1,:)
    ponts(cHull(t,2)+1,:)
    ponts(cHull(t,3)+1,:)
    ponts(cHull(t,1)+1,:)
    ]
    //id=color("green")
    plot2d(firstTri(:,1),firstTri(:,2))

    if (loud) then
        point1x = mean(firstTri(1:3,1)) +0.2*(mean(firstTri(1:2,1)) - mean(firstTri(1:3,1)));
        point1y = mean(firstTri(1:3,2)) +0.2*(mean(firstTri(1:2,2)) - mean(firstTri(1:3,2)));

        point2x = mean(firstTri(1:3,1)) +0.2*(mean(firstTri(2:3,1)) - mean(firstTri(1:3,1)));
        point2y = mean(firstTri(1:3,2)) +0.2*(mean(firstTri(2:3,2)) - mean(firstTri(1:3,2)));

        point3x = mean(firstTri(1:3,1)) +0.2*(mean(firstTri(3:4,1)) - mean(firstTri(1:3,1)));
        point3y = mean(firstTri(1:3,2)) +0.2*(mean(firstTri(3:4,2)) - mean(firstTri(1:3,2)));
        xstring(point1x,point1y,string(cHull(t,3+1)));
        xstring(point2x,point2y,string(cHull(t,3+2)));
        xstring(point3x,point3y,string(cHull(t,3+3)));

        //dx=-0.02
        //dy=-0.02
        dx=0;
        dy=0;
        xstring(mean(firstTri(1:3,1))+dx,mean(firstTri(1:3,2))+dy,string(t-1))

        id=color("blue")
        t=get("current_entity")
        t.font_foreground=id
    end
end
scatter(ponts(:,1),ponts(:,2),,"black",".")
dx=0.005
dy=-0.01
if (loud) then
    xstring(ponts(:,1)+dx,ponts(:,2)+dy,string(ponts(:,3)-1))
end

a=gca();
//a.axes_visible = ["off" "off" "off"];
a.box = "on"
a.data_bounds = [-0.1,-0.1;1.1,1.1]
a.data_bounds = [-5.1,-5.1;5.1,5.1]







scf(3)
clf(3)
ponts = zeros(154,3);
ponts(1,1)= 1.000; ponts(1,2)=-1.000; ponts(1,3)=1; ponts(1,4)= 0.000;
ponts(2,1)= 1.000; ponts(2,2)= 1.000; ponts(2,3)=2; ponts(2,4)= 0.000;
ponts(3,1)=-1.000; ponts(3,2)= 1.000; ponts(3,3)=3; ponts(3,4)= 0.000;
ponts(4,1)=-1.000; ponts(4,2)=-1.000; ponts(4,3)=4; ponts(4,4)= 0.000;
ponts(5,1)= 0.725; ponts(5,2)= 2.738; ponts(5,3)=5; ponts(5,4)= 0.000;
ponts(6,1)= 0.711; ponts(6,2)= 4.471; ponts(6,3)=6; ponts(6,4)= 0.000;
ponts(7,1)=-0.941; ponts(7,2)= 4.507; ponts(7,3)=7; ponts(7,4)= 0.000;
ponts(8,1)=-0.949; ponts(8,2)= 2.758; ponts(8,3)=8; ponts(8,4)= 0.000;
ponts(9,1)=-2.470; ponts(9,2)= 0.741; ponts(9,3)=9; ponts(9,4)= 0.000;
ponts(10,1)=-2.470; ponts(10,2)= 2.741; ponts(10,3)=10; ponts(10,4)= 0.000;
ponts(11,1)=-4.470; ponts(11,2)= 2.741; ponts(11,3)=11; ponts(11,4)= 0.933;
ponts(12,1)=-4.470; ponts(12,2)= 0.741; ponts(12,3)=12; ponts(12,4)= 0.933;
ponts(13,1)= 1.000; ponts(13,2)= 0.000; ponts(13,3)=13; ponts(13,4)= 0.000;
ponts(14,1)= 0.000; ponts(14,2)= 1.000; ponts(14,3)=14; ponts(14,4)= 0.000;
ponts(15,1)=-1.000; ponts(15,2)= 0.000; ponts(15,3)=15; ponts(15,4)= 0.000;
ponts(16,1)= 0.000; ponts(16,2)=-1.000; ponts(16,3)=16; ponts(16,4)= 0.000;
ponts(17,1)= 1.000; ponts(17,2)=-0.500; ponts(17,3)=17; ponts(17,4)= 0.000;
ponts(18,1)= 0.250; ponts(18,2)=-0.250; ponts(18,3)=18; ponts(18,4)= 0.000;
ponts(19,1)= 0.500; ponts(19,2)=-1.000; ponts(19,3)=19; ponts(19,4)= 0.000;
ponts(20,1)=-0.500; ponts(20,2)=-0.500; ponts(20,3)=20; ponts(20,4)= 0.000;
ponts(21,1)=-1.000; ponts(21,2)=-0.500; ponts(21,3)=21; ponts(21,4)= 0.000;
ponts(22,1)=-0.500; ponts(22,2)=-1.000; ponts(22,3)=22; ponts(22,4)= 0.000;
ponts(23,1)=-0.312; ponts(23,2)= 0.188; ponts(23,3)=23; ponts(23,4)= 0.000;
ponts(24,1)=-0.750; ponts(24,2)=-0.750; ponts(24,3)=24; ponts(24,4)= 0.000;
ponts(25,1)=-0.750; ponts(25,2)=-1.000; ponts(25,3)=25; ponts(25,4)= 0.000;
ponts(26,1)=-1.000; ponts(26,2)=-0.750; ponts(26,3)=26; ponts(26,4)= 0.000;
ponts(27,1)=-0.625; ponts(27,2)=-0.875; ponts(27,3)=27; ponts(27,4)= 0.000;
ponts(28,1)=-0.625; ponts(28,2)=-1.000; ponts(28,3)=28; ponts(28,4)= 0.000;
ponts(29,1)=-0.375; ponts(29,2)=-0.750; ponts(29,3)=29; ponts(29,4)= 0.000;
ponts(30,1)=-0.098; ponts(30,2)=-0.455; ponts(30,3)=30; ponts(30,4)= 0.000;
ponts(31,1)= 0.259; ponts(31,2)=-0.664; ponts(31,3)=31; ponts(31,4)= 0.000;
ponts(32,1)=-0.562; ponts(32,2)=-0.688; ponts(32,3)=32; ponts(32,4)= 0.000;
ponts(33,1)=-0.719; ponts(33,2)=-0.531; ponts(33,3)=33; ponts(33,4)= 0.000;
ponts(34,1)=-0.660; ponts(34,2)=-0.160; ponts(34,3)=34; ponts(34,4)= 0.000;
ponts(35,1)=-1.000; ponts(35,2)= 0.500; ponts(35,3)=35; ponts(35,4)= 0.000;
ponts(36,1)=-0.690; ponts(36,2)= 0.217; ponts(36,3)=36; ponts(36,4)= 0.000;
ponts(37,1)=-0.500; ponts(37,2)= 1.000; ponts(37,3)=37; ponts(37,4)= 0.000;
ponts(38,1)=-0.472; ponts(38,2)= 0.579; ponts(38,3)=38; ponts(38,4)= 0.000;
ponts(39,1)= 0.625; ponts(39,2)=-0.375; ponts(39,3)=39; ponts(39,4)= 0.000;
ponts(40,1)= 0.562; ponts(40,2)= 0.062; ponts(40,3)=40; ponts(40,4)= 0.000;
ponts(41,1)= 1.000; ponts(41,2)= 0.500; ponts(41,3)=41; ponts(41,4)= 0.000;
ponts(42,1)= 0.500; ponts(42,2)= 1.000; ponts(42,3)=42; ponts(42,4)= 0.000;
ponts(43,1)= 0.250; ponts(43,2)= 0.513; ponts(43,3)=43; ponts(43,4)= 0.000;
ponts(44,1)= 0.750; ponts(44,2)=-0.750; ponts(44,3)=44; ponts(44,4)= 0.000;
ponts(45,1)= 0.750; ponts(45,2)=-1.000; ponts(45,3)=45; ponts(45,4)= 0.000;
ponts(46,1)= 1.000; ponts(46,2)=-0.750; ponts(46,3)=46; ponts(46,4)= 0.000;
ponts(47,1)= 0.380; ponts(47,2)=-0.832; ponts(47,3)=47; ponts(47,4)= 0.000;
ponts(48,1)= 0.250; ponts(48,2)=-1.000; ponts(48,3)=48; ponts(48,4)= 0.000;
ponts(49,1)= 0.513; ponts(49,2)=-0.609; ponts(49,3)=49; ponts(49,4)= 0.000;
ponts(50,1)= 0.688; ponts(50,2)=-0.562; ponts(50,3)=50; ponts(50,4)= 0.000;
ponts(51,1)= 1.000; ponts(51,2)=-0.250; ponts(51,3)=51; ponts(51,4)= 0.000;
ponts(52,1)= 0.820; ponts(52,2)=-0.414; ponts(52,3)=52; ponts(52,4)= 0.000;
ponts(53,1)= 0.759; ponts(53,2)=-0.125; ponts(53,3)=53; ponts(53,4)= 0.000;
ponts(54,1)= 0.875; ponts(54,2)=-0.875; ponts(54,3)=54; ponts(54,4)= 0.000;
ponts(55,1)= 1.000; ponts(55,2)=-0.875; ponts(55,3)=55; ponts(55,4)= 0.000;
ponts(56,1)= 0.875; ponts(56,2)=-1.000; ponts(56,3)=56; ponts(56,4)= 0.000;
ponts(57,1)= 0.438; ponts(57,2)=-0.312; ponts(57,3)=57; ponts(57,4)= 0.000;
ponts(58,1)= 0.406; ponts(58,2)=-0.094; ponts(58,3)=58; ponts(58,4)= 0.000;
ponts(59,1)= 0.273; ponts(59,2)= 0.195; ponts(59,3)=59; ponts(59,4)= 0.000;
ponts(60,1)= 0.285; ponts(60,2)=-0.456; ponts(60,3)=60; ponts(60,4)= 0.000;
ponts(61,1)= 0.179; ponts(61,2)=-0.023; ponts(61,3)=61; ponts(61,4)= 0.000;
ponts(62,1)=-0.020; ponts(62,2)= 0.192; ponts(62,3)=62; ponts(62,4)= 0.000;
ponts(63,1)=-0.161; ponts(63,2)=-0.119; ponts(63,3)=63; ponts(63,4)= 0.000;
ponts(64,1)=-0.171; ponts(64,2)= 0.473; ponts(64,3)=64; ponts(64,4)= 0.000;
ponts(65,1)= 0.598; ponts(65,2)=-0.219; ponts(65,3)=65; ponts(65,4)= 0.000;
ponts(66,1)=-1.000; ponts(66,2)=-0.250; ponts(66,3)=66; ponts(66,4)= 0.000;
ponts(67,1)=-0.859; ponts(67,2)=-0.391; ponts(67,3)=67; ponts(67,4)= 0.000;
ponts(68,1)=-0.689; ponts(68,2)=-0.346; ponts(68,3)=68; ponts(68,4)= 0.000;
ponts(69,1)=-0.481; ponts(69,2)=-0.283; ponts(69,3)=69; ponts(69,4)= 0.000;
ponts(70,1)=-0.440; ponts(70,2)=-0.032; ponts(70,3)=70; ponts(70,4)= 0.000;
ponts(71,1)=-0.049; ponts(71,2)=-0.728; ponts(71,3)=71; ponts(71,4)= 0.000;
ponts(72,1)= 0.100; ponts(72,2)=-0.864; ponts(72,3)=72; ponts(72,4)= 0.000;
ponts(73,1)=-0.830; ponts(73,2)=-0.205; ponts(73,3)=73; ponts(73,4)= 0.000;
ponts(74,1)=-0.793; ponts(74,2)=-0.001; ponts(74,3)=74; ponts(74,4)= 0.000;
ponts(75,1)=-1.000; ponts(75,2)= 0.250; ponts(75,3)=75; ponts(75,4)= 0.000;
ponts(76,1)=-0.875; ponts(76,2)=-0.625; ponts(76,3)=76; ponts(76,4)= 0.000;
ponts(77,1)=-1.000; ponts(77,2)=-0.625; ponts(77,3)=77; ponts(77,4)= 0.000;
ponts(78,1)=-0.617; ponts(78,2)=-0.017; ponts(78,3)=78; ponts(78,4)= 0.000;
ponts(79,1)=-0.514; ponts(79,2)= 0.144; ponts(79,3)=79; ponts(79,4)= 0.000;
ponts(80,1)=-0.455; ponts(80,2)= 0.358; ponts(80,3)=80; ponts(80,4)= 0.000;
ponts(81,1)=-0.671; ponts(81,2)= 0.452; ponts(81,3)=81; ponts(81,4)= 0.000;
ponts(82,1)=-0.750; ponts(82,2)= 1.000; ponts(82,3)=82; ponts(82,4)= 0.000;
ponts(83,1)=-1.000; ponts(83,2)= 0.750; ponts(83,3)=83; ponts(83,4)= 0.000;
ponts(84,1)=-0.705; ponts(84,2)= 0.727; ponts(84,3)=84; ponts(84,4)= 0.000;
ponts(85,1)=-0.875; ponts(85,2)=-0.875; ponts(85,3)=85; ponts(85,4)= 0.000;
ponts(86,1)=-0.875; ponts(86,2)=-1.000; ponts(86,3)=86; ponts(86,4)= 0.000;
ponts(87,1)=-1.000; ponts(87,2)=-0.875; ponts(87,3)=87; ponts(87,4)= 0.000;
ponts(88,1)=-0.250; ponts(88,2)=-1.000; ponts(88,3)=88; ponts(88,4)= 0.000;
ponts(89,1)=-0.150; ponts(89,2)=-0.864; ponts(89,3)=89; ponts(89,4)= 0.000;
ponts(90,1)=-0.212; ponts(90,2)=-0.739; ponts(90,3)=90; ponts(90,4)= 0.000;
ponts(91,1)=-0.306; ponts(91,2)=-0.559; ponts(91,3)=91; ponts(91,4)= 0.000;
ponts(92,1)=-0.155; ponts(92,2)=-0.597; ponts(92,3)=92; ponts(92,4)= 0.000;
ponts(93,1)= 0.002; ponts(93,2)=-0.578; ponts(93,3)=93; ponts(93,4)= 0.000;
ponts(94,1)= 0.094; ponts(94,2)=-0.401; ponts(94,3)=94; ponts(94,4)= 0.000;
ponts(95,1)= 0.105; ponts(95,2)=-0.696; ponts(95,3)=95; ponts(95,4)= 0.000;
ponts(96,1)= 0.009; ponts(96,2)=-0.071; ponts(96,3)=96; ponts(96,4)= 0.000;
ponts(97,1)=-0.030; ponts(97,2)=-0.257; ponts(97,3)=97; ponts(97,4)= 0.000;
ponts(98,1)=-0.217; ponts(98,2)=-0.304; ponts(98,3)=98; ponts(98,4)= 0.000;
ponts(99,1)=-0.394; ponts(99,2)=-0.421; ponts(99,3)=99; ponts(99,4)= 0.000;
ponts(100,1)=-0.835; ponts(100,2)= 0.351; ponts(100,3)=100; ponts(100,4)= 0.000;
ponts(101,1)=-0.845; ponts(101,2)= 0.234; ponts(101,3)=101; ponts(101,4)= 0.000;
ponts(102,1)= 0.781; ponts(102,2)= 0.281; ponts(102,3)=102; ponts(102,4)= 0.000;
ponts(103,1)= 1.000; ponts(103,2)= 0.250; ponts(103,3)=103; ponts(103,4)= 0.000;
ponts(104,1)= 0.641; ponts(104,2)= 0.641; ponts(104,3)=104; ponts(104,4)= 0.000;
ponts(105,1)= 0.262; ponts(105,2)= 0.354; ponts(105,3)=105; ponts(105,4)= 0.000;
ponts(106,1)= 0.487; ponts(106,2)= 0.450; ponts(106,3)=106; ponts(106,4)= 0.000;
ponts(107,1)= 0.047; ponts(107,2)= 0.418; ponts(107,3)=107; ponts(107,4)= 0.000;
ponts(108,1)= 0.034; ponts(108,2)= 0.710; ponts(108,3)=108; ponts(108,4)= 0.000;
ponts(109,1)= 0.126; ponts(109,2)= 0.264; ponts(109,3)=109; ponts(109,4)= 0.000;
ponts(110,1)= 0.424; ponts(110,2)= 0.286; ponts(110,3)=110; ponts(110,4)= 0.000;
ponts(111,1)= 0.603; ponts(111,2)= 0.312; ponts(111,3)=111; ponts(111,4)= 0.000;
ponts(112,1)=-0.250; ponts(112,2)= 1.000; ponts(112,3)=112; ponts(112,4)= 0.000;
ponts(113,1)= 0.721; ponts(113,2)= 0.465; ponts(113,3)=113; ponts(113,4)= 0.000;
ponts(114,1)= 1.000; ponts(114,2)= 0.750; ponts(114,3)=114; ponts(114,4)= 0.000;
ponts(115,1)= 0.841; ponts(115,2)= 0.626; ponts(115,3)=115; ponts(115,4)= 0.000;
ponts(116,1)= 0.418; ponts(116,2)= 0.129; ponts(116,3)=116; ponts(116,4)= 0.000;
ponts(117,1)= 0.750; ponts(117,2)= 1.000; ponts(117,3)=117; ponts(117,4)= 0.000;
ponts(118,1)=-0.486; ponts(118,2)= 0.789; ponts(118,3)=118; ponts(118,4)= 0.000;
ponts(119,1)=-0.260; ponts(119,2)= 0.699; ponts(119,3)=119; ponts(119,4)= 0.000;
ponts(120,1)=-0.321; ponts(120,2)= 0.526; ponts(120,3)=120; ponts(120,4)= 0.000;
ponts(121,1)=-0.293; ponts(121,2)= 0.366; ponts(121,3)=121; ponts(121,4)= 0.000;
ponts(122,1)=-0.136; ponts(122,2)= 0.311; ponts(122,3)=122; ponts(122,4)= 0.000;
ponts(123,1)= 0.875; ponts(123,2)=-0.625; ponts(123,3)=123; ponts(123,4)= 0.000;
ponts(124,1)= 1.000; ponts(124,2)=-0.625; ponts(124,3)=124; ponts(124,4)= 0.000;
ponts(125,1)=-0.068; ponts(125,2)= 0.591; ponts(125,3)=125; ponts(125,4)= 0.000;
ponts(126,1)= 0.770; ponts(126,2)= 0.078; ponts(126,3)=126; ponts(126,4)= 0.000;
ponts(127,1)= 0.375; ponts(127,2)= 0.756; ponts(127,3)=127; ponts(127,4)= 0.000;
ponts(128,1)= 0.431; ponts(128,2)= 0.603; ponts(128,3)=128; ponts(128,4)= 0.000;
ponts(129,1)= 0.625; ponts(129,2)=-0.875; ponts(129,3)=129; ponts(129,4)= 0.000;
ponts(130,1)= 0.625; ponts(130,2)=-1.000; ponts(130,3)=130; ponts(130,4)= 0.000;
ponts(131,1)= 0.475; ponts(131,2)=-0.461; ponts(131,3)=131; ponts(131,4)= 0.000;
ponts(132,1)=-0.237; ponts(132,2)= 0.034; ponts(132,3)=132; ponts(132,4)= 0.000;
ponts(133,1)= 0.891; ponts(133,2)= 0.391; ponts(133,3)=133; ponts(133,4)= 0.000;
ponts(134,1)= 1.000; ponts(134,2)= 0.375; ponts(134,3)=134; ponts(134,4)= 0.000;
ponts(135,1)= 0.176; ponts(135,2)=-0.532; ponts(135,3)=135; ponts(135,4)= 0.000;
ponts(136,1)= 0.570; ponts(136,2)= 0.820; ponts(136,3)=136; ponts(136,4)= 0.000;
ponts(137,1)= 0.706; ponts(137,2)= 0.723; ponts(137,3)=137; ponts(137,4)= 0.000;
ponts(138,1)= 0.875; ponts(138,2)= 0.875; ponts(138,3)=138; ponts(138,4)= 0.000;
ponts(139,1)= 1.000; ponts(139,2)= 0.875; ponts(139,3)=139; ponts(139,4)= 0.000;
ponts(140,1)= 0.875; ponts(140,2)= 1.000; ponts(140,3)=140; ponts(140,4)= 0.000;
ponts(141,1)=-0.368; ponts(141,2)= 0.895; ponts(141,3)=141; ponts(141,4)= 0.000;
ponts(142,1)=-0.375; ponts(142,2)= 1.000; ponts(142,3)=142; ponts(142,4)= 0.000;
ponts(143,1)=-0.875; ponts(143,2)= 0.875; ponts(143,3)=143; ponts(143,4)= 0.000;
ponts(144,1)=-0.875; ponts(144,2)= 1.000; ponts(144,3)=144; ponts(144,4)= 0.000;
ponts(145,1)=-1.000; ponts(145,2)= 0.875; ponts(145,3)=145; ponts(145,4)= 0.000;
ponts(146,1)=-0.853; ponts(146,2)= 0.613; ponts(146,3)=146; ponts(146,4)= 0.000;
ponts(147,1)= 0.130; ponts(147,2)=-0.161; ponts(147,3)=147; ponts(147,4)= 0.000;
ponts(148,1)=-0.329; ponts(148,2)=-0.168; ponts(148,3)=148; ponts(148,4)= 0.000;
ponts(149,1)= 0.250; ponts(149,2)= 1.000; ponts(149,3)=149; ponts(149,4)= 0.000;
ponts(150,1)= 0.142; ponts(150,2)= 0.855; ponts(150,3)=150; ponts(150,4)= 0.000;
ponts(151,1)= 0.205; ponts(151,2)= 0.733; ponts(151,3)=151; ponts(151,4)= 0.000;
ponts(152,1)=-0.618; ponts(152,2)= 0.895; ponts(152,3)=152; ponts(152,4)= 0.000;
ponts(153,1)=-0.625; ponts(153,2)= 1.000; ponts(153,3)=153; ponts(153,4)= 0.000;
ponts(154,1)=-0.130; ponts(154,2)= 0.849; ponts(154,3)=154; ponts(154,4)= 0.000;
cHull = [ ]
cHull=[cHull; 16 50 51 -1 7 42 ]
cHull=[cHull; 45 53 54 52 67 -1 ]
cHull=[cHull; 38 51 52 61 7 86 ]
cHull=[cHull; 52 39 64 153 28 86 ]
cHull=[cHull; 32 19 67 16 93 91 ]
cHull=[cHull; 30 48 59 56 70 26 ]
cHull=[cHull; 46 30 71 56 138 99 ]
cHull=[cHull; 51 50 52 0 60 2 ]
cHull=[cHull; 17 57 60 11 176 229 ]
cHull=[cHull; 66 67 72 91 87 101 ]
cHull=[cHull; 32 66 75 91 19 106 ]
cHull=[cHull; 17 56 57 69 41 8 ]
cHull=[cHull; 23 31 32 13 16 106 ]
cHull=[cHull; 23 26 31 20 29 12 ]
cHull=[cHull; 69 22 78 205 94 111 ]
cHull=[cHull; 33 68 69 92 146 109 ]
cHull=[cHull; 31 19 32 30 4 12 ]
cHull=[cHull; 24 84 85 18 122 -1 ]
cHull=[cHull; 24 23 84 20 120 17 ]
cHull=[cHull; 66 20 75 88 107 10 ]
cHull=[cHull; 23 24 26 18 22 13 ]
cHull=[cHull; 21 26 27 23 22 -1 ]
cHull=[cHull; 26 24 27 20 -1 21 ]
cHull=[cHull; 26 21 28 21 124 29 ]
cHull=[cHull; 29 90 91 78 96 55 ]
cHull=[cHull; 19 28 90 30 31 130 ]
cHull=[cHull; 30 59 134 5 75 211 ]
cHull=[cHull; 18 46 47 57 99 -1 ]
cHull=[cHull; 39 57 64 177 41 3 ]
cHull=[cHull; 26 28 31 23 30 13 ]
cHull=[cHull; 28 19 31 25 16 29 ]
cHull=[cHull; 28 89 90 127 96 25 ]
cHull=[cHull; 29 92 93 55 133 48 ]
cHull=[cHull; 68 19 98 93 130 206 ]
cHull=[cHull; 73 35 100 90 37 150 ]
cHull=[cHull; 35 80 99 36 227 37 ]
cHull=[cHull; 35 79 80 47 40 35 ]
cHull=[cHull; 35 99 100 35 95 34 ]
cHull=[cHull; 39 110 115 63 54 177 ]
cHull=[cHull; 37 79 119 40 112 184 ]
cHull=[cHull; 79 37 80 39 119 36 ]
cHull=[cHull; 57 56 64 11 85 28 ]
cHull=[cHull; 16 51 122 0 59 190 ]
cHull=[cHull; 105 109 110 166 54 168 ]
cHull=[cHull; 107 149 153 236 170 242 ]
cHull=[cHull; 101 110 125 169 63 195 ]
cHull=[cHull; 105 103 127 152 213 199 ]
cHull=[cHull; 35 78 79 110 94 36 ]
cHull=[cHull; 29 93 96 32 136 81 ]
cHull=[cHull; 38 56 130 85 64 204 ]
cHull=[cHull; 28 87 88 124 97 127 ]
cHull=[cHull; 43 48 128 58 53 201 ]
cHull=[cHull; 45 43 53 189 65 1 ]
cHull=[cHull; 48 46 128 56 57 51 ]
cHull=[cHull; 110 109 115 43 72 38 ]
cHull=[cHull; 29 91 92 24 76 32 ]
cHull=[cHull; 30 46 48 6 53 5 ]
cHull=[cHull; 46 18 128 27 202 53 ]
cHull=[cHull; 48 43 49 51 62 203 ]
cHull=[cHull; 51 49 122 61 62 42 ]
cHull=[cHull; 50 12 52 -1 159 7 ]
cHull=[cHull; 38 49 51 204 59 2 ]
cHull=[cHull; 49 43 122 58 189 59 ]
cHull=[cHull; 110 39 125 38 153 45 ]
cHull=[cHull; 56 59 130 69 70 49 ]
cHull=[cHull; 43 44 53 201 68 52 ]
cHull=[cHull; 0 53 55 67 68 -1 ]
cHull=[cHull; 53 0 54 66 -1 1 ]
cHull=[cHull; 53 44 55 65 -1 66 ]
cHull=[cHull; 56 17 59 11 134 64 ]
cHull=[cHull; 59 48 130 5 203 64 ]
cHull=[cHull; 62 95 131 142 82 143 ]
cHull=[cHull; 109 58 115 167 158 54 ]
cHull=[cHull; 40 114 132 174 171 208 ]
cHull=[cHull; 104 106 108 161 151 165 ]
cHull=[cHull; 59 93 134 134 133 26 ]
cHull=[cHull; 91 70 92 131 132 55 ]
cHull=[cHull; 106 42 107 161 237 193 ]
cHull=[cHull; 90 29 97 24 81 145 ]
cHull=[cHull; 63 119 120 183 112 83 ]
cHull=[cHull; 1 137 138 219 218 -1 ]
cHull=[cHull; 29 96 97 48 129 78 ]
cHull=[cHull; 95 61 131 139 140 71 ]
cHull=[cHull; 63 120 121 79 185 188 ]
cHull=[cHull; 116 137 139 154 219 -1 ]
cHull=[cHull; 56 38 64 49 86 41 ]
cHull=[cHull; 38 52 64 2 3 85 ]
cHull=[cHull; 67 33 72 92 103 9 ]
cHull=[cHull; 65 20 66 -1 19 101 ]
cHull=[cHull; 25 75 76 105 107 -1 ]
cHull=[cHull; 35 73 77 34 108 110 ]
cHull=[cHull; 66 32 67 10 4 9 ]
cHull=[cHull; 33 67 68 87 93 15 ]
cHull=[cHull; 67 19 68 4 33 92 ]
cHull=[cHull; 78 22 79 14 113 47 ]
cHull=[cHull; 99 74 100 148 149 37 ]
cHull=[cHull; 90 89 91 31 131 24 ]
cHull=[cHull; 87 15 88 -1 125 50 ]
cHull=[cHull; 15 47 71 -1 99 125 ]
cHull=[cHull; 47 46 71 27 6 98 ]
cHull=[cHull; 14 65 72 -1 101 102 ]
cHull=[cHull; 65 66 72 88 9 100 ]
cHull=[cHull; 14 72 73 100 103 150 ]
cHull=[cHull; 72 33 73 87 108 102 ]
cHull=[cHull; 118 63 124 183 192 241 ]
cHull=[cHull; 25 23 75 120 106 89 ]
cHull=[cHull; 23 32 75 12 10 105 ]
cHull=[cHull; 75 20 76 19 -1 89 ]
cHull=[cHull; 73 33 77 103 109 90 ]
cHull=[cHull; 33 69 77 15 111 108 ]
cHull=[cHull; 35 77 78 90 111 47 ]
cHull=[cHull; 77 69 78 109 14 110 ]
cHull=[cHull; 119 79 120 39 113 79 ]
cHull=[cHull; 79 22 120 94 185 112 ]
cHull=[cHull; 117 118 140 182 235 221 ]
cHull=[cHull; 2 142 143 226 225 -1 ]
cHull=[cHull; 83 142 145 163 118 228 ]
cHull=[cHull; 82 142 144 118 226 -1 ]
cHull=[cHull; 142 82 145 117 147 116 ]
cHull=[cHull; 80 37 83 40 179 228 ]
cHull=[cHull; 23 25 84 105 123 18 ]
cHull=[cHull; 3 84 86 122 123 -1 ]
cHull=[cHull; 84 3 85 121 -1 17 ]
cHull=[cHull; 84 25 86 120 -1 121 ]
cHull=[cHull; 28 21 87 23 -1 50 ]
cHull=[cHull; 15 71 88 98 126 97 ]
cHull=[cHull; 71 70 88 137 128 125 ]
cHull=[cHull; 28 88 89 50 128 31 ]
cHull=[cHull; 88 70 89 126 131 127 ]
cHull=[cHull; 96 62 97 142 232 81 ]
cHull=[cHull; 19 90 98 25 145 33 ]
cHull=[cHull; 89 70 91 128 76 96 ]
cHull=[cHull; 92 70 94 76 137 210 ]
cHull=[cHull; 93 92 134 32 210 75 ]
cHull=[cHull; 59 17 93 69 141 75 ]
cHull=[cHull; 95 96 146 142 136 230 ]
cHull=[cHull; 96 93 146 48 141 135 ]
cHull=[cHull; 70 71 94 126 138 132 ]
cHull=[cHull; 71 30 94 6 211 137 ]
cHull=[cHull; 60 61 95 160 82 230 ]
cHull=[cHull; 61 22 131 186 205 82 ]
cHull=[cHull; 93 17 146 134 229 136 ]
cHull=[cHull; 95 62 96 71 129 135 ]
cHull=[cHull; 62 131 147 71 144 232 ]
cHull=[cHull; 131 69 147 205 146 143 ]
cHull=[cHull; 90 97 98 78 231 130 ]
cHull=[cHull; 69 68 147 15 206 144 ]
cHull=[cHull; 82 34 145 -1 223 118 ]
cHull=[cHull; 34 74 99 -1 95 223 ]
cHull=[cHull; 74 14 100 -1 150 95 ]
cHull=[cHull; 14 73 100 102 34 149 ]
cHull=[cHull; 106 61 108 187 160 74 ]
cHull=[cHull; 103 105 112 46 168 172 ]
cHull=[cHull; 39 52 125 3 159 63 ]
cHull=[cHull; 116 136 137 178 216 84 ]
cHull=[cHull; 116 41 135 -1 175 178 ]
cHull=[cHull; 42 127 150 199 197 237 ]
cHull=[cHull; 42 104 105 161 166 199 ]
cHull=[cHull; 58 60 115 164 176 72 ]
cHull=[cHull; 52 12 125 60 194 153 ]
cHull=[cHull; 61 60 108 139 164 151 ]
cHull=[cHull; 104 42 106 157 77 74 ]
cHull=[cHull; 81 142 151 225 163 239 ]
cHull=[cHull; 142 83 151 116 224 162 ]
cHull=[cHull; 60 58 108 158 165 160 ]
cHull=[cHull; 58 104 108 167 74 164 ]
cHull=[cHull; 105 104 109 157 167 43 ]
cHull=[cHull; 104 58 109 165 72 166 ]
cHull=[cHull; 105 110 112 43 169 152 ]
cHull=[cHull; 110 101 112 45 173 168 ]
cHull=[cHull; 149 13 153 234 180 44 ]
cHull=[cHull; 114 112 132 172 173 73 ]
cHull=[cHull; 103 112 114 152 171 215 ]
cHull=[cHull; 112 101 132 169 207 171 ]
cHull=[cHull; 40 113 114 -1 217 73 ]
cHull=[cHull; 41 126 135 233 212 155 ]
cHull=[cHull; 60 57 115 8 177 158 ]
cHull=[cHull; 57 39 115 28 38 176 ]
cHull=[cHull; 116 135 136 155 214 154 ]
cHull=[cHull; 83 37 117 119 182 224 ]
cHull=[cHull; 13 111 153 -1 181 170 ]
cHull=[cHull; 111 140 153 222 235 180 ]
cHull=[cHull; 117 37 118 179 184 114 ]
cHull=[cHull; 63 118 119 104 184 79 ]
cHull=[cHull; 118 37 119 182 39 183 ]
cHull=[cHull; 120 22 121 113 186 83 ]
cHull=[cHull; 22 61 121 140 187 185 ]
cHull=[cHull; 61 106 121 151 188 186 ]
cHull=[cHull; 106 63 121 192 83 187 ]
cHull=[cHull; 43 45 122 52 191 62 ]
cHull=[cHull; 16 122 123 42 191 -1 ]
cHull=[cHull; 122 45 123 189 -1 190 ]
cHull=[cHull; 63 106 124 188 193 104 ]
cHull=[cHull; 106 107 124 77 242 192 ]
cHull=[cHull; 12 102 125 -1 195 159 ]
cHull=[cHull; 102 101 125 207 45 194 ]
cHull=[cHull; 126 148 149 233 234 198 ]
cHull=[cHull; 127 126 150 212 198 156 ]
cHull=[cHull; 126 149 150 196 236 197 ]
cHull=[cHull; 42 105 127 157 46 156 ]
cHull=[cHull; 44 128 129 201 202 -1 ]
cHull=[cHull; 44 43 128 65 51 200 ]
cHull=[cHull; 128 18 129 57 -1 200 ]
cHull=[cHull; 48 49 130 58 204 70 ]
cHull=[cHull; 49 38 130 61 49 203 ]
cHull=[cHull; 22 69 131 14 144 140 ]
cHull=[cHull; 68 98 147 33 231 146 ]
cHull=[cHull; 101 102 132 195 209 173 ]
cHull=[cHull; 40 132 133 73 209 -1 ]
cHull=[cHull; 132 102 133 207 -1 208 ]
cHull=[cHull; 92 94 134 132 211 133 ]
cHull=[cHull; 94 30 134 138 26 210 ]
cHull=[cHull; 126 127 135 197 213 175 ]
cHull=[cHull; 127 103 135 46 214 212 ]
cHull=[cHull; 135 103 136 213 215 178 ]
cHull=[cHull; 103 114 136 172 216 214 ]
cHull=[cHull; 136 114 137 215 217 154 ]
cHull=[cHull; 114 113 137 174 218 216 ]
cHull=[cHull; 137 113 138 217 -1 80 ]
cHull=[cHull; 137 1 139 80 -1 84 ]
cHull=[cHull; 36 140 141 221 222 -1 ]
cHull=[cHull; 36 117 140 238 114 220 ]
cHull=[cHull; 140 111 141 181 -1 220 ]
cHull=[cHull; 34 99 145 148 227 147 ]
cHull=[cHull; 83 117 151 179 238 163 ]
cHull=[cHull; 142 81 143 162 -1 115 ]
cHull=[cHull; 142 2 144 115 -1 117 ]
cHull=[cHull; 99 80 145 35 228 223 ]
cHull=[cHull; 80 83 145 119 116 227 ]
cHull=[cHull; 17 60 146 8 230 141 ]
cHull=[cHull; 60 95 146 139 135 229 ]
cHull=[cHull; 98 97 147 145 232 206 ]
cHull=[cHull; 97 62 147 129 143 231 ]
cHull=[cHull; 126 41 148 175 -1 196 ]
cHull=[cHull; 148 13 149 -1 170 196 ]
cHull=[cHull; 140 118 153 114 241 181 ]
cHull=[cHull; 149 107 150 44 237 198 ]
cHull=[cHull; 107 42 150 77 156 236 ]
cHull=[cHull; 117 36 151 221 240 224 ]
cHull=[cHull; 81 151 152 162 240 -1 ]
cHull=[cHull; 151 36 152 238 -1 239 ]
cHull=[cHull; 118 124 153 104 242 235 ]
cHull=[cHull; 124 107 153 193 44 241 ]

loud = %T
for t=1:(size(cHull)(1))
    firstTri=[
    ponts(cHull(t,1)+1,:)
    ponts(cHull(t,2)+1,:)
    ponts(cHull(t,3)+1,:)
    ponts(cHull(t,1)+1,:)
    ]
    //id=color("green")
    plot2d(firstTri(:,1),firstTri(:,2))

    if (loud) then
        point1x = mean(firstTri(1:3,1)) +0.2*(mean(firstTri(1:2,1)) - mean(firstTri(1:3,1)));
        point1y = mean(firstTri(1:3,2)) +0.2*(mean(firstTri(1:2,2)) - mean(firstTri(1:3,2)));

        point2x = mean(firstTri(1:3,1)) +0.2*(mean(firstTri(2:3,1)) - mean(firstTri(1:3,1)));
        point2y = mean(firstTri(1:3,2)) +0.2*(mean(firstTri(2:3,2)) - mean(firstTri(1:3,2)));

        point3x = mean(firstTri(1:3,1)) +0.2*(mean(firstTri(3:4,1)) - mean(firstTri(1:3,1)));
        point3y = mean(firstTri(1:3,2)) +0.2*(mean(firstTri(3:4,2)) - mean(firstTri(1:3,2)));
        xstring(point1x,point1y,string(cHull(t,3+1)));
        xstring(point2x,point2y,string(cHull(t,3+2)));
        xstring(point3x,point3y,string(cHull(t,3+3)));

        //dx=-0.02
        //dy=-0.02
        dx=0;
        dy=0;
        xstring(mean(firstTri(1:3,1))+dx,mean(firstTri(1:3,2))+dy,string(t-1))

        id=color("blue")
        t=get("current_entity")
        t.font_foreground=id
    end
end
scatter(ponts(:,1),ponts(:,2),,"black",".")
dx=0.005
dy=-0.01
if (loud) then
    xstring(ponts(:,1)+dx,ponts(:,2)+dy,string(ponts(:,3)-1))
end

a=gca();
//a.axes_visible = ["off" "off" "off"];
a.box = "on"
//a.data_bounds = [-0.1,-0.1;1.1,1.1]
a.data_bounds = [-5.1,-5.1;5.1,5.1]





//N=10;               //Number of points
//points=rand(N,2);

if (%F) then
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

//scatter(points(:,1),points(:,2),,"red",".")
