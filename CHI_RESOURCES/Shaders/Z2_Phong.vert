#version 330
//#########################################################
//#                                                       #
//#                  PHONG SHADER (VERTEX)                #
//#                                                       #
//#########################################################
//=============================================== Incoming vertex data
layout (location = 0) in      vec3 vVertex;
layout (location = 1) in      vec3 vNormal;
layout (location = 2) in      vec3 vTexCoord0;
layout (location = 3) in      vec3 vTangent;
layout (location = 4) in      vec3 vBinorm;

//=============================================== Interpolations
out     vec3        texCoord0;
out     vec3        normal;
out     vec3        bi_norm;
out     vec3        tangent;

out     vec3        pixelWorldSpace;     //Absolute world space
out     vec3        pixelViewSpace;      //Pixel to eye

//=============================================== Uniforms
uniform mat4        projectionMatrix;   // 0Projection matrix
uniform mat4        viewMatrix;         // 1View matrix
uniform mat4        modelMatrix;        // 2Model matrix
uniform mat4        normalMatrix;       // 3Normal matrix
uniform mat4        textureMatrix;      // 4Texture matrix

uniform int         lightCount;         // 5Number of lights
uniform int         lightEnabled[10];   // 6 0=false, 1=true
uniform int         lightType[10];      // 7 0=Sun,1=point,2=spot
uniform vec3        lightPosition[10];  // 8Doubles as Sun-direction
uniform vec4        lightParams[10];    // 9Misc. light params
uniform ivec4       lightBParams[10];   //10Light bool-like params
uniform vec4        lightColor[10];     //11Light color

uniform vec4        materialDiffuse;    //12Material diffuse color
uniform vec4        materialAmbient;    //13Material ambient color
uniform vec4        materialSpecular;   //14Material specular color
uniform float       materialShininess;  //15Material shininess

uniform sampler2D   textureIndex[6];    //16Texture indices
uniform int         textureFlags[12];   //17Texture enabled flags

uniform ivec4       vobjectID;          //18Object ID

//=============================================== Main function
void main()
{
        texCoord0   = (textureMatrix*vec4(vTexCoord0.x,vTexCoord0.y,vTexCoord0.z,1.0)).xyz;
        normal      = normalize(normalMatrix*vec4(vNormal.x,vNormal.y,vNormal.z,1.0)).xyz;
        tangent     = normalize(normalMatrix*vec4(vTangent.x,vTangent.y,vTangent.z,1.0)).xyz;
        bi_norm     = normalize(normalMatrix*vec4(vBinorm.x,vBinorm.y,vBinorm.z,1.0)).xyz;

        pixelWorldSpace = (modelMatrix*vec4(vVertex.x,vVertex.y,vVertex.z,1.0)).xyz;
        pixelViewSpace  = (viewMatrix*modelMatrix*vec4(vVertex.x,vVertex.y,vVertex.z,1.0)).xyz;

        gl_Position =projectionMatrix*viewMatrix*modelMatrix*vec4(vVertex.x,vVertex.y,vVertex.z,1.0);
}