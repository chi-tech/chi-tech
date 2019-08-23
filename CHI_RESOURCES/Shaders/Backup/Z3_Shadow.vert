#version 130
//#########################################################
//#                                                       #
//#         SHADOWMAP GENERATION SHADER (VERTEX)          #
//#                                                       #
//#########################################################
in      vec4 vVertex;
varying vec3 shadowCoord;

uniform mat4 lWorldMatrix;

varying     vec4 fragPos;
varying     vec3 lightPos;


void main()
{        
        vec3 cameraSpaceVertex = (gl_ModelViewMatrix*vec4(vVertex.x,vVertex.y,vVertex.z,1.0)).xyz;
        lightPos=(gl_LightSource[0].position).xyz;
        fragPos=vec4(cameraSpaceVertex,1.0);
        //
        //shadowCoord = (gl_ModelViewMatrix*vVertex).xyz - gl_LightSource[0].position.xyz;
        shadowCoord = (lWorldMatrix*vVertex).xyz - gl_LightSource[0].position.xyz;
        gl_Position = gl_ModelViewProjectionMatrix*vVertex;
}