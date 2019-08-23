#version 330
//#########################################################
//#                                                       #
//#                  FLAT SHADER (FRAGMENT)               #
//#                                                       #
//#########################################################
//=============================================== Incoming vertex data



in  vec4 pixelColor;

//=============================================== Uniforms
uniform mat4        projectionMatrix;   // 0Projection matrix
uniform mat4        viewMatrix;         // 1View matrix
uniform mat4        modelMatrix;        // 2Model matrix

//=============================================== Color output definitions
layout(location = 0) out vec4 color0;


//############################################### Main function
void main (void)
{
    color0 = pixelColor;
}
//###############################################

