#version 330
//#########################################################
//#                                                       #
//#                PICKING SHADER (VERTEX)                #
//#                                                       #
//#########################################################
//=============================================== Incoming vertex data
layout (location = 0) in      vec3 vVertex;




//=============================================== Uniforms
uniform mat4        projectionMatrix;   // 0Projection matrix
uniform mat4        viewMatrix;         // 1View matrix
uniform mat4        modelMatrix;        // 2Model matrix

uniform ivec4       vobjectID;          //18Object ID

//=============================================== Main function
void main()
{
        gl_Position =projectionMatrix*viewMatrix*modelMatrix*vec4(vVertex.x,vVertex.y,vVertex.z,1.0);
}