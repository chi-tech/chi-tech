#version 330
//#########################################################
//#                                                       #
//#                PICKING SHADER (FRAGMENT)              #
//#                                                       #
//#########################################################
//=============================================== Incoming vertex data






//=============================================== Uniforms
uniform mat4        projectionMatrix;   // 0Projection matrix
uniform mat4        viewMatrix;         // 1View matrix
uniform mat4        modelMatrix;        // 2Model matrix

uniform ivec4       vobjectID;          //18Object ID

//=============================================== Color output definitions
layout(location = 0) out vec4 color0;


//############################################### Main function
void main (void)
{
    color0.r = float(vobjectID.x/255.0);
    color0.g = float(vobjectID.y/255.0);
    color0.b = float(vobjectID.z/255.0);
    color0.a = 1.0;
}
//###############################################

