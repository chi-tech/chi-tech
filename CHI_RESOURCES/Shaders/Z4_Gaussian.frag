#version 130
//#########################################################
//#                                                       #
//#           GAUSSIAN BLUR SHADER (FRAGMENT)             #
//#                                                       #
//#########################################################

uniform         float       bSigma;
uniform         float       bSize;
uniform         sampler2D   blurSample;

const           float       pi = 3.14159265f;
const           float       pixelPerSide = 3.0f;
const           vec2        multiHoriVector = vec2(1.0f, 0.0f);
//const           vec2        multiVertVector = vec2(0.0f, 1.0f);

void main() {
    vec3 iGaussian;

    iGaussian.x = 1.0f/(sqrt(2.0f * pi)*bSigma);
    iGaussian.y = exp(-0.5f/(bSigma*bSigma));
    iGaussian.z = iGaussian.y*iGaussian.x;

    vec4 avgValue = vec4(0.0f, 0.0f, 0.0f, 0.0f);
    float iSum = 0.0f;

    avgValue += texture2D(blurSample, gl_TexCoord[0].xy)*iGaussian.x;
    iSum += iGaussian.x;
    iGaussian.xy *= iGaussian.xy;

    for (int i = 1; i <= pixelPerSide; i++)
    {
      avgValue += texture2D(blurSample, gl_TexCoord[0].xy - i * bSize * multiHoriVector) * iGaussian.x;
      avgValue += texture2D(blurSample, gl_TexCoord[0].xy + i * bSize * multiHoriVector) * iGaussian.x;
      iSum += 2 * iGaussian.x;
      iGaussian.xy *= iGaussian.yz;
    }

    gl_FragColor = avgValue / iSum;
}
