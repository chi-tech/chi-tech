#version 130
//#########################################################
//#                                                       #
//#         SELECTION MAP GENERATION SHADER (FRAGMENT)    #
//#                                                       #
//#########################################################
varying     vec4 objectID;

void main (void)
{
        gl_FragColor.r = float(objectID.x);
        gl_FragColor.g = float(objectID.y);
        gl_FragColor.b = float(objectID.z);
        gl_FragColor.a = 1.0;
        //gl_FragColor=vec4(1.0,0.0,0.0,1.0);
}

