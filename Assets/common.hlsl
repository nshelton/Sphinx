
float2 _Resolution;
int _Level;
float4x4 _CameraToWorld;
float _Threshold;

const float PI = 3.14159265359;


struct Ray 
{
    float3 origin;
    float3 direction;
};

Ray getRay(float2 uv) {
    
    uv -= 0.5;
    uv *= float2(1, _Resolution.y/_Resolution.x);
    float3 dir = normalize(float3(uv, 1)) ;

    Ray result;

    result.origin = float3(_CameraToWorld._m03, _CameraToWorld._m13, _CameraToWorld._m23);
    result.direction = mul(_CameraToWorld, float4(dir, 0)).xyz;

    return result;
}

float getPixelSize(float d) {
    return pow(2, _Level) * d * _Threshold;
}

float rand(float2 co)
{
    return frac(sin(dot(co.xy, float2(12.9898, 78.233))) * 43758.5453) - 0.5;
}

#define SEED 0
// PRNG function
float nrand(float2 uv, float salt)
{
    uv += float2(salt, SEED);
    return frac(sin(dot(uv, float2(12.9898, 78.233))) * 43758.5453);
}

//Random number [0:1] without sine
#define HASHSCALE1 .1031
float hash(float p)
{
	float3 p3  = frac((float3)p * HASHSCALE1);
    p3 += dot(p3, p3.yzx + 19.19);
    return frac((p3.x + p3.y) * p3.z);
}

float3 randomSphereDir(float2 uv)
{
    float u = nrand(uv, 10) * 2 - 1;
    float theta = nrand(uv, 11) * PI * 2;
    float u2 = sqrt(1 - u * u);
    return float3(u2 * cos(theta), u2 * sin(theta), u);
}

float3 randomHemisphereDir(float3 dir, float i)
{
	float3 v = randomSphereDir( float2(rand(i+1.), rand(i+2.)) );
	return v * sign(dot(v, dir));
}
 
