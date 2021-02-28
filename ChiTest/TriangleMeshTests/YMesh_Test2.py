import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import mpl_toolkits.mplot3d as a3
from matplotlib import colors
import scipy as sp
xyz=np.zeros((99,3))
xyz[0][0]=-20.000000;  xyz[0][1]=-20.000000;  xyz[0][2]=0.000000;
xyz[1][0]=20.000000;  xyz[1][1]=-20.000000;  xyz[1][2]=0.000000;
xyz[2][0]=-20.000000;  xyz[2][1]=20.000000;  xyz[2][2]=0.000000;
xyz[3][0]=20.000000;  xyz[3][1]=20.000000;  xyz[3][2]=0.000000;
xyz[4][0]=-10.567614;  xyz[4][1]=-10.341662;  xyz[4][2]=0.000004;
xyz[5][0]=-10.706662;  xyz[5][1]=12.253562;  xyz[5][2]=0.000008;
xyz[6][0]=-4.520560;  xyz[6][1]=12.392611;  xyz[6][2]=0.000008;
xyz[7][0]=-4.172941;  xyz[7][1]=-10.411185;  xyz[7][2]=0.000008;
xyz[8][0]=-4.346750;  xyz[8][1]=8.585411;  xyz[8][2]=0.000008;
xyz[9][0]=-4.259846;  xyz[9][1]=-7.279806;  xyz[9][2]=0.000008;
xyz[10][0]=10.011423;  xyz[10][1]=-7.491186;  xyz[10][2]=0.000008;
xyz[11][0]=9.733327;  xyz[11][1]=8.568805;  xyz[11][2]=0.000008;
xyz[12][0]=-10.637138;  xyz[12][1]=8.531524;  xyz[12][2]=0.000006;
xyz[13][0]=-10.602377;  xyz[13][1]=-7.402124;  xyz[13][2]=0.000005;
xyz[14][0]=-20.000000;  xyz[14][1]=-7.435355;  xyz[14][2]=0.000000;
xyz[15][0]=-20.000000;  xyz[15][1]=-10.273048;  xyz[15][2]=0.000000;
xyz[16][0]=-20.000000;  xyz[16][1]=8.474646;  xyz[16][2]=0.000000;
xyz[17][0]=-20.000000;  xyz[17][1]=12.417519;  xyz[17][2]=0.000000;
xyz[18][0]=-4.595316;  xyz[18][1]=20.000000;  xyz[18][2]=0.000000;
xyz[19][0]=-10.695234;  xyz[19][1]=20.000000;  xyz[19][2]=0.000000;
xyz[20][0]=10.000000;  xyz[20][1]=20.000000;  xyz[20][2]=0.000000;
xyz[21][0]=9.866663;  xyz[21][1]=12.372454;  xyz[21][2]=0.000004;
xyz[22][0]=9.733326;  xyz[22][1]=-10.758802;  xyz[22][2]=0.000008;
xyz[23][0]=-10.706663;  xyz[23][1]=-19.796896;  xyz[23][2]=0.000004;
xyz[24][0]=-3.893331;  xyz[24][1]=-19.866419;  xyz[24][2]=0.000004;
xyz[25][0]=10.219995;  xyz[25][1]=-19.727369;  xyz[25][2]=0.000004;
xyz[26][0]=20.161896;  xyz[26][1]=-10.341662;  xyz[26][2]=0.000000;
xyz[27][0]=19.883802;  xyz[27][1]=-7.421662;  xyz[27][2]=0.000000;
xyz[28][0]=19.814276;  xyz[28][1]=8.499281;  xyz[28][2]=0.000004;
xyz[29][0]=19.466656;  xyz[29][1]=12.462135;  xyz[29][2]=0.000004;
xyz[30][0]=-20.000000;  xyz[30][1]=-15.136524;  xyz[30][2]=0.000000;
xyz[31][0]=-15.353332;  xyz[31][1]=-19.898449;  xyz[31][2]=0.000002;
xyz[32][0]=20.080948;  xyz[32][1]=-15.170832;  xyz[32][2]=0.000000;
xyz[33][0]=-15.347616;  xyz[33][1]=20.000000;  xyz[33][2]=0.000000;
xyz[34][0]=-10.671900;  xyz[34][1]=10.392543;  xyz[34][2]=0.000007;
xyz[35][0]=-7.613611;  xyz[35][1]=12.323086;  xyz[35][2]=0.000008;
xyz[36][0]=-4.433655;  xyz[36][1]=10.489011;  xyz[36][2]=0.000008;
xyz[37][0]=-7.370277;  xyz[37][1]=-10.376424;  xyz[37][2]=0.000006;
xyz[38][0]=-4.216393;  xyz[38][1]=-8.845495;  xyz[38][2]=0.000008;
xyz[39][0]=-4.303298;  xyz[39][1]=0.652803;  xyz[39][2]=0.000008;
xyz[40][0]=2.875789;  xyz[40][1]=-7.385496;  xyz[40][2]=0.000008;
xyz[41][0]=9.872375;  xyz[41][1]=0.538809;  xyz[41][2]=0.000008;
xyz[42][0]=2.693288;  xyz[42][1]=8.577108;  xyz[42][2]=0.000008;
xyz[43][0]=-10.619757;  xyz[43][1]=0.564700;  xyz[43][2]=0.000005;
xyz[44][0]=-10.584995;  xyz[44][1]=-8.871893;  xyz[44][2]=0.000004;
xyz[45][0]=-7.491944;  xyz[45][1]=8.558468;  xyz[45][2]=0.000007;
xyz[46][0]=-7.431111;  xyz[46][1]=-7.340965;  xyz[46][2]=0.000006;
xyz[47][0]=-20.000000;  xyz[47][1]=0.519645;  xyz[47][2]=0.000000;
xyz[48][0]=-20.000000;  xyz[48][1]=-8.854202;  xyz[48][2]=0.000000;
xyz[49][0]=-20.000000;  xyz[49][1]=10.446082;  xyz[49][2]=0.000000;
xyz[50][0]=-15.283806;  xyz[50][1]=-10.307355;  xyz[50][2]=0.000002;
xyz[51][0]=-15.301188;  xyz[51][1]=-7.418740;  xyz[51][2]=0.000002;
xyz[52][0]=-15.318569;  xyz[52][1]=8.503085;  xyz[52][2]=0.000003;
xyz[53][0]=-20.000000;  xyz[53][1]=16.208759;  xyz[53][2]=0.000000;
xyz[54][0]=-15.353332;  xyz[54][1]=12.335541;  xyz[54][2]=0.000004;
xyz[55][0]=2.702342;  xyz[55][1]=20.000000;  xyz[55][2]=0.000000;
xyz[56][0]=-7.645276;  xyz[56][1]=20.000000;  xyz[56][2]=0.000000;
xyz[57][0]=15.000000;  xyz[57][1]=20.000000;  xyz[57][2]=0.000000;
xyz[58][0]=-10.700949;  xyz[58][1]=16.126781;  xyz[58][2]=0.000004;
xyz[59][0]=-4.557938;  xyz[59][1]=16.196306;  xyz[59][2]=0.000004;
xyz[60][0]=9.933331;  xyz[60][1]=16.186226;  xyz[60][2]=0.000002;
xyz[61][0]=9.799995;  xyz[61][1]=10.470629;  xyz[61][2]=0.000006;
xyz[62][0]=2.673052;  xyz[62][1]=12.382531;  xyz[62][2]=0.000006;
xyz[63][0]=2.780193;  xyz[63][1]=-10.584993;  xyz[63][2]=0.000008;
xyz[64][0]=9.872375;  xyz[64][1]=-9.124994;  xyz[64][2]=0.000008;
xyz[65][0]=-7.299997;  xyz[65][1]=-19.831657;  xyz[65][2]=0.000004;
xyz[66][0]=3.163332;  xyz[66][1]=-19.796894;  xyz[66][2]=0.000004;
xyz[67][0]=-10.637138;  xyz[67][1]=-15.069280;  xyz[67][2]=0.000004;
xyz[68][0]=-4.033136;  xyz[68][1]=-15.138803;  xyz[68][2]=0.000006;
xyz[69][0]=9.976661;  xyz[69][1]=-15.243087;  xyz[69][2]=0.000006;
xyz[70][0]=15.109997;  xyz[70][1]=-19.863684;  xyz[70][2]=0.000002;
xyz[71][0]=20.022850;  xyz[71][1]=-8.881662;  xyz[71][2]=0.000000;
xyz[72][0]=19.849039;  xyz[72][1]=0.538809;  xyz[72][2]=0.000002;
xyz[73][0]=19.640467;  xyz[73][1]=10.480708;  xyz[73][2]=0.000004;
xyz[74][0]=19.733330;  xyz[74][1]=16.231068;  xyz[74][2]=0.000002;
xyz[75][0]=14.666660;  xyz[75][1]=12.417295;  xyz[75][2]=0.000004;
xyz[76][0]=14.773800;  xyz[76][1]=8.534042;  xyz[76][2]=0.000006;
xyz[77][0]=14.947613;  xyz[77][1]=-7.456424;  xyz[77][2]=0.000004;
xyz[78][0]=14.947611;  xyz[78][1]=-10.550232;  xyz[78][2]=0.000004;
xyz[79][0]=-15.309878;  xyz[79][1]=0.542172;  xyz[79][2]=0.000003;
xyz[80][0]=2.784539;  xyz[80][1]=0.595806;  xyz[80][2]=0.000008;
xyz[81][0]=-7.461528;  xyz[81][1]=0.608751;  xyz[81][2]=0.000006;
xyz[82][0]=-7.629443;  xyz[82][1]=16.161543;  xyz[82][2]=0.000004;
xyz[83][0]=-15.350474;  xyz[83][1]=16.167770;  xyz[83][2]=0.000002;
xyz[84][0]=-15.318569;  xyz[84][1]=-15.102901;  xyz[84][2]=0.000002;
xyz[85][0]=-15.292497;  xyz[85][1]=-8.863048;  xyz[85][2]=0.000002;
xyz[86][0]=14.860706;  xyz[86][1]=0.538809;  xyz[86][2]=0.000005;
xyz[87][0]=-7.400694;  xyz[87][1]=-8.858694;  xyz[87][2]=0.000006;
xyz[88][0]=-15.335951;  xyz[88][1]=10.419312;  xyz[88][2]=0.000003;
xyz[89][0]=15.028803;  xyz[89][1]=-15.206958;  xyz[89][2]=0.000003;
xyz[90][0]=-7.335137;  xyz[90][1]=-15.104041;  xyz[90][2]=0.000005;
xyz[91][0]=14.720230;  xyz[91][1]=10.475669;  xyz[91][2]=0.000005;
xyz[92][0]=2.971762;  xyz[92][1]=-15.190944;  xyz[92][2]=0.000006;
xyz[93][0]=14.947613;  xyz[93][1]=-9.003328;  xyz[93][2]=0.000004;
xyz[94][0]=2.687697;  xyz[94][1]=16.191265;  xyz[94][2]=0.000003;
xyz[95][0]=14.833330;  xyz[95][1]=16.208647;  xyz[95][2]=0.000002;
xyz[96][0]=2.827991;  xyz[96][1]=-8.985245;  xyz[96][2]=0.000008;
xyz[97][0]=2.683170;  xyz[97][1]=10.479820;  xyz[97][2]=0.000007;
xyz[98][0]=-7.552777;  xyz[98][1]=10.440777;  xyz[98][2]=0.000007;
face_numverts=np.zeros((80,2))   
face_vertindi=np.zeros((80,10))
face_numverts[0,0]=4   
face_vertindi[0][0]=98
face_vertindi[0][1]=36
face_vertindi[0][2]=6
face_vertindi[0][3]=35
face_numverts[1,0]=4   
face_vertindi[1][0]=97
face_vertindi[1][1]=61
face_vertindi[1][2]=21
face_vertindi[1][3]=62
face_numverts[2,0]=4   
face_vertindi[2][0]=96
face_vertindi[2][1]=40
face_vertindi[2][2]=9
face_vertindi[2][3]=38
face_numverts[3,0]=4   
face_vertindi[3][0]=95
face_vertindi[3][1]=75
face_vertindi[3][2]=29
face_vertindi[3][3]=74
face_numverts[4,0]=4   
face_vertindi[4][0]=51
face_vertindi[4][1]=79
face_vertindi[4][2]=47
face_vertindi[4][3]=14
face_numverts[5,0]=4   
face_vertindi[5][0]=94
face_vertindi[5][1]=62
face_vertindi[5][2]=21
face_vertindi[5][3]=60
face_numverts[6,0]=4   
face_vertindi[6][0]=93
face_vertindi[6][1]=77
face_vertindi[6][2]=10
face_vertindi[6][3]=64
face_numverts[7,0]=4   
face_vertindi[7][0]=92
face_vertindi[7][1]=68
face_vertindi[7][2]=24
face_vertindi[7][3]=66
face_numverts[8,0]=4   
face_vertindi[8][0]=91
face_vertindi[8][1]=73
face_vertindi[8][2]=29
face_vertindi[8][3]=75
face_numverts[9,0]=4   
face_vertindi[9][0]=90
face_vertindi[9][1]=67
face_vertindi[9][2]=23
face_vertindi[9][3]=65
face_numverts[10,0]=4   
face_vertindi[10][0]=89
face_vertindi[10][1]=69
face_vertindi[10][2]=25
face_vertindi[10][3]=70
face_numverts[11,0]=4   
face_vertindi[11][0]=88
face_vertindi[11][1]=49
face_vertindi[11][2]=16
face_vertindi[11][3]=52
face_numverts[12,0]=4   
face_vertindi[12][0]=87
face_vertindi[12][1]=38
face_vertindi[12][2]=9
face_vertindi[12][3]=46
face_numverts[13,0]=4   
face_vertindi[13][0]=86
face_vertindi[13][1]=72
face_vertindi[13][2]=28
face_vertindi[13][3]=76
face_numverts[14,0]=4   
face_vertindi[14][0]=85
face_vertindi[14][1]=48
face_vertindi[14][2]=15
face_vertindi[14][3]=50
face_numverts[15,0]=4   
face_vertindi[15][0]=84
face_vertindi[15][1]=31
face_vertindi[15][2]=23
face_vertindi[15][3]=67
face_numverts[16,0]=4   
face_vertindi[16][0]=83
face_vertindi[16][1]=53
face_vertindi[16][2]=17
face_vertindi[16][3]=54
face_numverts[17,0]=4   
face_vertindi[17][0]=82
face_vertindi[17][1]=56
face_vertindi[17][2]=19
face_vertindi[17][3]=58
face_numverts[18,0]=4   
face_vertindi[18][0]=81
face_vertindi[18][1]=46
face_vertindi[18][2]=9
face_vertindi[18][3]=39
face_numverts[19,0]=4   
face_vertindi[19][0]=80
face_vertindi[19][1]=41
face_vertindi[19][2]=11
face_vertindi[19][3]=42
face_numverts[20,0]=4   
face_vertindi[20][0]=58
face_vertindi[20][1]=83
face_vertindi[20][2]=54
face_vertindi[20][3]=5
face_numverts[21,0]=4   
face_vertindi[21][0]=43
face_vertindi[21][1]=13
face_vertindi[21][2]=46
face_vertindi[21][3]=81
face_numverts[22,0]=4   
face_vertindi[22][0]=9
face_vertindi[22][1]=40
face_vertindi[22][2]=80
face_vertindi[22][3]=39
face_numverts[23,0]=4   
face_vertindi[23][0]=45
face_vertindi[23][1]=81
face_vertindi[23][2]=39
face_vertindi[23][3]=8
face_numverts[24,0]=4   
face_vertindi[24][0]=13
face_vertindi[24][1]=43
face_vertindi[24][2]=79
face_vertindi[24][3]=51
face_numverts[25,0]=4   
face_vertindi[25][0]=50
face_vertindi[25][1]=84
face_vertindi[25][2]=67
face_vertindi[25][3]=4
face_numverts[26,0]=4   
face_vertindi[26][0]=43
face_vertindi[26][1]=12
face_vertindi[26][2]=52
face_vertindi[26][3]=79
face_numverts[27,0]=4   
face_vertindi[27][0]=40
face_vertindi[27][1]=10
face_vertindi[27][2]=41
face_vertindi[27][3]=80
face_numverts[28,0]=4   
face_vertindi[28][0]=39
face_vertindi[28][1]=80
face_vertindi[28][2]=42
face_vertindi[28][3]=8
face_numverts[29,0]=4   
face_vertindi[29][0]=35
face_vertindi[29][1]=82
face_vertindi[29][2]=58
face_vertindi[29][3]=5
face_numverts[30,0]=4   
face_vertindi[30][0]=12
face_vertindi[30][1]=43
face_vertindi[30][2]=81
face_vertindi[30][3]=45
face_numverts[31,0]=4   
face_vertindi[31][0]=33
face_vertindi[31][1]=2
face_vertindi[31][2]=53
face_vertindi[31][3]=83
face_numverts[32,0]=4   
face_vertindi[32][0]=59
face_vertindi[32][1]=18
face_vertindi[32][2]=56
face_vertindi[32][3]=82
face_numverts[33,0]=4   
face_vertindi[33][0]=15
face_vertindi[33][1]=30
face_vertindi[33][2]=84
face_vertindi[33][3]=50
face_numverts[34,0]=4   
face_vertindi[34][0]=13
face_vertindi[34][1]=51
face_vertindi[34][2]=85
face_vertindi[34][3]=44
face_numverts[35,0]=4   
face_vertindi[35][0]=19
face_vertindi[35][1]=33
face_vertindi[35][2]=83
face_vertindi[35][3]=58
face_numverts[36,0]=4   
face_vertindi[36][0]=51
face_vertindi[36][1]=14
face_vertindi[36][2]=48
face_vertindi[36][3]=85
face_numverts[37,0]=4   
face_vertindi[37][0]=30
face_vertindi[37][1]=0
face_vertindi[37][2]=31
face_vertindi[37][3]=84
face_numverts[38,0]=4   
face_vertindi[38][0]=44
face_vertindi[38][1]=85
face_vertindi[38][2]=50
face_vertindi[38][3]=4
face_numverts[39,0]=4   
face_vertindi[39][0]=6
face_vertindi[39][1]=59
face_vertindi[39][2]=82
face_vertindi[39][3]=35
face_numverts[40,0]=4   
face_vertindi[40][0]=79
face_vertindi[40][1]=52
face_vertindi[40][2]=16
face_vertindi[40][3]=47
face_numverts[41,0]=4   
face_vertindi[41][0]=41
face_vertindi[41][1]=86
face_vertindi[41][2]=76
face_vertindi[41][3]=11
face_numverts[42,0]=4   
face_vertindi[42][0]=10
face_vertindi[42][1]=77
face_vertindi[42][2]=86
face_vertindi[42][3]=41
face_numverts[43,0]=4   
face_vertindi[43][0]=77
face_vertindi[43][1]=27
face_vertindi[43][2]=72
face_vertindi[43][3]=86
face_numverts[44,0]=4   
face_vertindi[44][0]=44
face_vertindi[44][1]=87
face_vertindi[44][2]=46
face_vertindi[44][3]=13
face_numverts[45,0]=4   
face_vertindi[45][0]=4
face_vertindi[45][1]=37
face_vertindi[45][2]=87
face_vertindi[45][3]=44
face_numverts[46,0]=4   
face_vertindi[46][0]=37
face_vertindi[46][1]=7
face_vertindi[46][2]=38
face_vertindi[46][3]=87
face_numverts[47,0]=4   
face_vertindi[47][0]=34
face_vertindi[47][1]=88
face_vertindi[47][2]=52
face_vertindi[47][3]=12
face_numverts[48,0]=4   
face_vertindi[48][0]=5
face_vertindi[48][1]=54
face_vertindi[48][2]=88
face_vertindi[48][3]=34
face_numverts[49,0]=4   
face_vertindi[49][0]=54
face_vertindi[49][1]=17
face_vertindi[49][2]=49
face_vertindi[49][3]=88
face_numverts[50,0]=4   
face_vertindi[50][0]=32
face_vertindi[50][1]=89
face_vertindi[50][2]=70
face_vertindi[50][3]=1
face_numverts[51,0]=4   
face_vertindi[51][0]=26
face_vertindi[51][1]=78
face_vertindi[51][2]=89
face_vertindi[51][3]=32
face_numverts[52,0]=4   
face_vertindi[52][0]=78
face_vertindi[52][1]=22
face_vertindi[52][2]=69
face_vertindi[52][3]=89
face_numverts[53,0]=4   
face_vertindi[53][0]=68
face_vertindi[53][1]=90
face_vertindi[53][2]=65
face_vertindi[53][3]=24
face_numverts[54,0]=4   
face_vertindi[54][0]=7
face_vertindi[54][1]=37
face_vertindi[54][2]=90
face_vertindi[54][3]=68
face_numverts[55,0]=4   
face_vertindi[55][0]=37
face_vertindi[55][1]=4
face_vertindi[55][2]=67
face_vertindi[55][3]=90
face_numverts[56,0]=4   
face_vertindi[56][0]=61
face_vertindi[56][1]=91
face_vertindi[56][2]=75
face_vertindi[56][3]=21
face_numverts[57,0]=4   
face_vertindi[57][0]=11
face_vertindi[57][1]=76
face_vertindi[57][2]=91
face_vertindi[57][3]=61
face_numverts[58,0]=4   
face_vertindi[58][0]=76
face_vertindi[58][1]=28
face_vertindi[58][2]=73
face_vertindi[58][3]=91
face_numverts[59,0]=4   
face_vertindi[59][0]=69
face_vertindi[59][1]=92
face_vertindi[59][2]=66
face_vertindi[59][3]=25
face_numverts[60,0]=4   
face_vertindi[60][0]=22
face_vertindi[60][1]=63
face_vertindi[60][2]=92
face_vertindi[60][3]=69
face_numverts[61,0]=4   
face_vertindi[61][0]=63
face_vertindi[61][1]=7
face_vertindi[61][2]=68
face_vertindi[61][3]=92
face_numverts[62,0]=4   
face_vertindi[62][0]=78
face_vertindi[62][1]=93
face_vertindi[62][2]=64
face_vertindi[62][3]=22
face_numverts[63,0]=4   
face_vertindi[63][0]=26
face_vertindi[63][1]=71
face_vertindi[63][2]=93
face_vertindi[63][3]=78
face_numverts[64,0]=4   
face_vertindi[64][0]=71
face_vertindi[64][1]=27
face_vertindi[64][2]=77
face_vertindi[64][3]=93
face_numverts[65,0]=4   
face_vertindi[65][0]=55
face_vertindi[65][1]=94
face_vertindi[65][2]=60
face_vertindi[65][3]=20
face_numverts[66,0]=4   
face_vertindi[66][0]=18
face_vertindi[66][1]=59
face_vertindi[66][2]=94
face_vertindi[66][3]=55
face_numverts[67,0]=4   
face_vertindi[67][0]=59
face_vertindi[67][1]=6
face_vertindi[67][2]=62
face_vertindi[67][3]=94
face_numverts[68,0]=4   
face_vertindi[68][0]=57
face_vertindi[68][1]=95
face_vertindi[68][2]=74
face_vertindi[68][3]=3
face_numverts[69,0]=4   
face_vertindi[69][0]=20
face_vertindi[69][1]=60
face_vertindi[69][2]=95
face_vertindi[69][3]=57
face_numverts[70,0]=4   
face_vertindi[70][0]=60
face_vertindi[70][1]=21
face_vertindi[70][2]=75
face_vertindi[70][3]=95
face_numverts[71,0]=4   
face_vertindi[71][0]=63
face_vertindi[71][1]=96
face_vertindi[71][2]=38
face_vertindi[71][3]=7
face_numverts[72,0]=4   
face_vertindi[72][0]=22
face_vertindi[72][1]=64
face_vertindi[72][2]=96
face_vertindi[72][3]=63
face_numverts[73,0]=4   
face_vertindi[73][0]=64
face_vertindi[73][1]=10
face_vertindi[73][2]=40
face_vertindi[73][3]=96
face_numverts[74,0]=4   
face_vertindi[74][0]=36
face_vertindi[74][1]=97
face_vertindi[74][2]=62
face_vertindi[74][3]=6
face_numverts[75,0]=4   
face_vertindi[75][0]=8
face_vertindi[75][1]=42
face_vertindi[75][2]=97
face_vertindi[75][3]=36
face_numverts[76,0]=4   
face_vertindi[76][0]=42
face_vertindi[76][1]=11
face_vertindi[76][2]=61
face_vertindi[76][3]=97
face_numverts[77,0]=4   
face_vertindi[77][0]=34
face_vertindi[77][1]=98
face_vertindi[77][2]=35
face_vertindi[77][3]=5
face_numverts[78,0]=4   
face_vertindi[78][0]=12
face_vertindi[78][1]=45
face_vertindi[78][2]=98
face_vertindi[78][3]=34
face_numverts[79,0]=4   
face_vertindi[79][0]=45
face_vertindi[79][1]=8
face_vertindi[79][2]=36
face_vertindi[79][3]=98
fig = plt.figure(1)
ax = a3.Axes3D(plt.figure(1))
num_faces=80

for f in range(0,num_faces):
    vcount = int(face_numverts[f,0])
    verts=np.zeros((vcount+1,3))
    
    for v in range(0,vcount):
        index = int(face_vertindi[f][v])
        verts[v][0] = xyz[index][0]
        verts[v][1] = xyz[index][1]
        verts[v][2] = xyz[index][2]
    
    index = int(face_vertindi[f][0])
    verts[vcount][0] = xyz[index][0]
    verts[vcount][1] = xyz[index][1]
    verts[vcount][2] = xyz[index][2]
    
    #ax.plot(verts[:,0],verts[:,1],verts[:,2],'k-',linewidth=0.5)
    pol = a3.art3d.Poly3DCollection([verts],linewidths=0.5)
    pol.set_edgecolor('k')
    pol.set_facecolor(colors.rgb2hex([0.8,0.8,0.8]))
    if (int(face_numverts[f,1])==1):
        pol.set_facecolor(colors.rgb2hex([1.0,0.0,0.0]))
    ax.add_collection3d(pol)

xmax = np.max(xyz[:,0])
xmin = np.min(xyz[:,0])
ymax = np.max(xyz[:,1])
ymin = np.min(xyz[:,1])
zmax = np.max(xyz[:,2])
zmin = np.min(xyz[:,2])

ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])
ax.set_zlim([zmin, zmax])

ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.view_init(elev=90.0,azim=0.0)
plt.show()
