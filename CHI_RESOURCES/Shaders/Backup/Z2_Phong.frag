#version 130
//#########################################################
//#                                                       #
//#                  PHONG SHADER (FRAGMENT)              #
//#                                                       #
//#########################################################
in      vec3 normal;
in      vec3 lightDir[3];
in      vec3 eyeVec;
in      vec3 shadowCoord[3];
in      vec3 texCoord0;
in      vec3 tangent;
in      vec3 bi_norm;
in      vec4 objectID;
in      vec4 fragPos;
in      vec3 lightPos;

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
    float bias = 0.0;
    //================================================= Omnidirectional lights
    for (int k=0; k<1; k++)
    {
            lightDistance = length(fragPos.xyz-lightPos);//length(shadowCoord[k]);
            //lightDistance = length(fragPos.xyz-(gl_LightSource[k].position).xyz);

            L = normalize(lightDir[k]);
            //bias=0.005*tan(acos(dot(N,L)));
            //bias=0.2*tan(acos(dot(N,L)));
            bias=0.1;
            if ((lightDistance < lightCutOff[k]) && (lightFlags[k] == 1))
            {
                    if (lightShadowFlags[k] == 1)
                    {
                            //shadowDepth  = textureCube(shadowMapIndex[k],shadowCoord[k]).b;
                            vec4 shade=textureCube(shadowMapIndex[0],normalize(shadowCoord[k]));
                            shadowDepth  = 10.0*shade.x+10.0*shade.y+10.0*shade.z;
                            if (shadowDepth>=30.0) {shadowDepth=1000.0;}
                            //shadowDepth = shade.z;
                            //shadowDepth=shade.r;
                            //pixelDepth = clamp((lightDistance-20.0)/30.0,0.0,1.0);
                            //pixelDepth  =shade.z;//clamp(lightDistance/10.0,0.0,1.0);
                            pixelDepth = lightDistance;

                            if (pixelDepth <(shadowDepth+bias)) {shadowFactor = 1.0;}
                            else                                 {shadowFactor = 0.0;}


                    }
                    else {shadowFactor = 1.0;}

                    attenuation = 1.0/(gl_LightSource[k].constantAttenuation+
                                    gl_LightSource[k].linearAttenuation*lightDistance+
                                    gl_LightSource[k].quadraticAttenuation*lightDistance*lightDistance);


                    //L = normalize(lightDir[k]);
                    diffuseIntensity = clamp(max(0.0,dot(N,L)),0.0,1.0);
                    diffuseColor += shadowFactor*attenuation*diffuseIntensity*gl_LightSource[k].diffuse*materialDiffuse;
                    diffuseColor.a=materialDiffuse.a;
                  //diffuseColor=vec4(pixelDepth,pixelDepth,pixelDepth,1.0);
                    R = (reflect(-L,N));
                    specularIntensity = clamp(max(0.0,dot(E,R)),0.0,1.0);
                    specularFactor = pow(specularIntensity,gl_FrontMaterial.shininess);
                    specularColor += shadowFactor*attenuation*gl_FrontLightProduct[k].specular*specularFactor;

                    //diffuseColor = vec4(0.0,0.0,0.5,1.0);
            }
    }
    //vec4 outColor=ambientColor + diffuseColor + specularColor;
    gl_FragData[0] = clamp(ambientColor + diffuseColor + specularColor,0.0,1.0);
    //if (ambientColor.a<0.1)
    //{
    //    gl_FragData[0] = vec4(ambientColor.a,ambientColor.a,ambientColor.a,1.0);
    //}
    //gl_FragData[0] = ambientColor;
    //gl_FragData[0] = vec4(outColor.a*outColor.xyz+(1.0-outColor.a)*gl_FragData[0].xyz,1.0);
    //gl_FragData[0] = vec4(ambientColor.xyz,0.0)+ diffuseColor;

}
