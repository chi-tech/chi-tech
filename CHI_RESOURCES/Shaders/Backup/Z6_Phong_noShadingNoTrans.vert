#version 130
//#########################################################
//#                                                       #
//#                  PHONG SHADER (VERTEX)                #
//#                                                       #
//#########################################################
in      vec3 vVertex;
in      vec3 vNormal;
in      vec3 vTexCoord0;
in      vec3 vTangent;
in      vec3 vBinorm;

varying     vec3 normal;
varying     vec3 eyeVec;
varying     vec3 texCoord0;
varying     vec3 tangent;
varying     vec3 bi_norm;
varying     vec4 objectID;
varying     vec4 fragPos;

uniform mat4  cameraMatrix;
uniform mat4  lWorldMatrix;
uniform mat4  textureMatrix;
uniform ivec4 vobjectID;

void main()
{       
        texCoord0 = (textureMatrix*vec4(vTexCoord0.x,vTexCoord0.y,vTexCoord0.z,1.0)).xyz;
        normal  = normalize(gl_NormalMatrix*vNormal);
        tangent = normalize(gl_NormalMatrix*vTangent);
        bi_norm = normalize(gl_NormalMatrix*vBinorm);

        vec3 cameraSpaceVertex = (gl_ModelViewMatrix*vec4(vVertex.x,vVertex.y,vVertex.z,1.0)).xyz;
        vec3 worldSpaceVertex  = (lWorldMatrix*vec4(vVertex.x,vVertex.y,vVertex.z,1.0)).xyz;

        fragPos=vec4(worldSpaceVertex,1.0);

        eyeVec = -cameraSpaceVertex;
             
        gl_Position =gl_ModelViewProjectionMatrix*vec4(vVertex.x,vVertex.y,vVertex.z,1.0);
        
        objectID.r = float(vobjectID.x)/255.0;
        objectID.g = float(vobjectID.y)/255.0;
        objectID.b = float(vobjectID.z)/255.0;
        objectID.a = float(vobjectID.a);
}