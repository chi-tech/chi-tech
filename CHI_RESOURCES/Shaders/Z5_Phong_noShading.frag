#version 130
//#########################################################
//#                                                       #
//#                  PHONG SHADER (FRAGMENT)              #
//#                                                       #
//#########################################################
in      vec3 normal;
in      vec3 eyeVec;
in      vec3 texCoord0;
in      vec3 tangent;
in      vec3 bi_norm;
in      vec4 objectID;
in      vec4 fragPos;

uniform samplerCube shadowMapIndex[3];
uniform sampler2D   textureIndex[6];
uniform int         textureFlags[6];
uniform int         lightShadowFlags[3];
uniform int         lightFlags[3];
uniform float       lightCutOff[3];



mat3 tangentMatrix;


void main (void)
{
    
    vec3 N = normalize(normal); //Normal vector
    vec3 L = vec3(0.0,0.0,0.0); //Light direction 
    vec3 R = vec3(0.0,0.0,0.0); //Reflection vector
    vec3 E = normalize(eyeVec); //Eye vector
    
    tangentMatrix    = mat3(normalize(tangent), normalize(bi_norm), N);
    
    vec4 materialDiffuse = gl_FrontMaterial.diffuse;
    
    
    
    vec4 ambientColor       = gl_FrontMaterial.ambient; //gl_LightSource[0].ambient*gl_FrontMaterial.ambient;
    vec4 diffuseColor       = vec4(0.0,0.0,0.0,0.0);
    vec4 specularColor      = vec4(0.0,0.0,0.0,0.0);
    float shadowDepth       = 1.0;
    float diffuseIntensity  = 0.0;
    float specularIntensity = 0.0;
    float specularFactor    = 0.0;
    float attenuation       = 1.0;
    float lightDistance     = 0.0;
    float pixelDepth        = 0.0;
    float shadowFactor      = 1.0;
    
    //================================================= Setting up textures
    if (textureFlags[0] == 1) {ambientColor    = texture2D(textureIndex[0],texCoord0.st);}
    if (textureFlags[1] == 1) {materialDiffuse = texture2D(textureIndex[1],texCoord0.st);}
    if (textureFlags[3] == 1) {N    = normalize(texture2D(textureIndex[3],texCoord0.st).xyz*2.0 -1.0);
                            N.y  = N.y*-1.0;
                            N    = normalize(tangentMatrix*N);}
    
    
    gl_FragData[0] = clamp(ambientColor + materialDiffuse,0.0,1.0);
    gl_FragData[0] = vec4(materialDiffuse.r,materialDiffuse.g,materialDiffuse.b,materialDiffuse.a);
}

