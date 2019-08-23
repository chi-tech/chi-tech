#version 130
//#########################################################
//#                                                       #
//#         SELECTION MAP GENERATION SHADER (VERTEX)      #
//#                                                       #
//#########################################################
in          vec4 vVertex;
uniform     ivec4 vobjectID;
varying     vec4 objectID;


void main()
{       
    objectID = vec4(float(vobjectID.x)/255.0,float(vobjectID.y)/255.0,float(vobjectID.z)/255.0,0.0);
    
     gl_Position = gl_ModelViewProjectionMatrix*vVertex;             
}