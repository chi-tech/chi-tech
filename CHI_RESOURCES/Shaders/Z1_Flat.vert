#version 330
//#########################################################
//#                                                       #
//#                  FLAT SHADER (VERTEX)                 #
//#                                                       #
//#########################################################
//=============================================== Incoming vertex data
layout (location = 0) in      vec3 vVertex;
layout (location = 1) in      vec4 vColor;

out vec4 pixelColor;

//=============================================== Uniforms
uniform mat4        projectionMatrix;   // 0Projection matrix
uniform mat4        viewMatrix;         // 1View matrix
uniform mat4        modelMatrix;        // 2Model matrix

//=============================================== Main function
void main()
{
    pixelColor = vColor;
        gl_Position =projectionMatrix*viewMatrix*modelMatrix*vec4(vVertex.x,vVertex.y,vVertex.z,1.0);
}
