#include "foundation.cginc"

float4 u_paramA;
float4 u_paramB;
float4 u_paramC;
float4 u_paramD;
float4 _Quaternion;
float _LEVELS;
float _FractalType;
  
float2 pseudo_knightyan(float3 p)
{
    float3 CSize = u_paramA.xyz;
    float DEfactor = 1.;
    float orbit = 0;

    for (int i = 0; i < 6; i++)
    {

        float3 start = p;
        p = 2. * clamp(p, -CSize, CSize) - p;
        p = rotateX(p, u_paramB.x);
        p = rotateY(p, u_paramB.y);
        p = rotateZ(p, u_paramB.z);
        float k = max(u_paramC.x / dot(p, p), 1.);
        p *= k;
        DEfactor *= k + 0.05;

        orbit += length(start - p);
    }

    float rxy = length(p.xy);
    float ds = max(rxy - 0.92784, abs(rxy * p.z) / length(p)) / DEfactor;

    return float2(ds, orbit);
}

float2 tglad_variant(float3 z0)
{
    // z0 = modc(z0, 2.0);
    float mr = u_paramD.x, mxr = u_paramD.y;

    float4 scale = (float4) (u_paramD.z), p0 = u_paramA.xyzz;
    float4 z = float4(z0, 1.0);
    float orbit = 0;
    for (int n = 0; n < _LEVELS; n++)
    {
        float3 start = z;
        z.xyz = clamp(z.xyz, -u_paramB.x, u_paramB.x) * 2.0 - z.xyz;
        z *= scale / clamp(dot(z.xyz, z.xyz), mr, mxr);
        z += p0;
        orbit += length(start - z);
    }
    float dS = (length(max(abs(z.xyz) - u_paramC.xyz, 0.0)) - 0) / z.w;
    return float2(dS, orbit);
}


float2 tglad(float3 z0)
{
    // z0 = modc(z0, 2.0);

    float mr = 0.25, mxr = 1.0;
    float4 scale = float4(-3.12, -3.12, -3.12, 3.12), p0 = u_paramA.xyzz;
    float4 z = float4(z0, 1.0);
    float orbit = 0;

    for (int n = 0; n < _LEVELS; n++)
    {
        float3 start = z.xyz;

        z.xyz = clamp(z.xyz, -u_paramB.x, u_paramB.x) * 2.0 - z.xyz;
        z *= scale / clamp(dot(z.xyz, z.xyz), mr, mxr);
        z += p0;
        orbit += length(start - z.xyz);
        

    }

    float dS = (length(max(abs(z.xyz) - float3(1.2, 49.0, 1.4), 0.0)) - 0.06) / z.w;
    return float2(dS, orbit);
}

void sphereFold(inout float3 z, inout float dz)
{

    float fixedRadius2 = u_paramA.x;
    float minRadius2 = u_paramA.y;

    float r2 = dot(z, z);
    if (r2 < minRadius2)
    {
        // linear inner scaling
        float temp = (fixedRadius2 / minRadius2);
        z *= temp;
        dz *= temp;
    }
    else if (r2 < fixedRadius2)
    {
        // this is the actual sphere inversion
        float temp = (fixedRadius2 / r2);
        z *= temp;
        dz *= temp;
    }
}

float4 fromtwovectors(float3 u, float3 v)
{
    u = normalize(u);
    v = normalize(v);
    float m = sqrt(2.f + 2.f * dot(u, v));
    float3 w = (1.f / m) * cross(u, v);
    return float4(w.x, w.y, w.z, 0.5f * m);
}



float2 hartverdrahtet(float3 p) {

    float3 cs = u_paramA.xyz;
    float fs = u_paramA.w;

    float3 fc = u_paramB.xyz;
    float fu = u_paramB.w;
    float orbit = 0.0;

    float  fd = 0.763;
   float dEfactor=1.;
   //int fractal_iterations = 12;
   for(int i=0;i<_LEVELS;i++){
    float3 start = p;
      //box folding
      p=2.*clamp(p, -cs, cs)-p;
      //inversion
      float k=max(fs/dot(p,p),1.);
      p*=k;
      dEfactor*=k;
      //julia seed
      p+=fc;

      orbit += length(start - p);
   }
   //call basic shape and scale its DE
   //need to adjust fractal_distancemult with non zero julia seed
   float rxy=length(p.xy)-fu;
   //distance from pos to the pseudo kleinian basic shape ...
   return float2(fd*max(rxy,abs(length(p.xy)*p.z)/sqrt(dot(p,p)))/abs(dEfactor), orbit);
}

/*
// distance function from Hartverdrahtet
// ( http://www.pouet.net/prod.php?which=59086 )
float2 hartverdrahtet(float3 f)
{
    float3 cs = u_paramA.xyz;
    float fs = u_paramA.w;

    float3 fc = u_paramB.xyz;
    float fu = u_paramB.w;
    float orbit = 0.0;
 
    float v = 1.;
    for (int i = 0; i < _LEVELS; i++)
    {
        float3 start = f;

        f = 2. * clamp(f, -cs, cs) - f;
        float c = max(fs / dot(f, f), 1.);
        //f = rotateX(f, u_paramC.x);
       // f = rotateY(f, u_paramC.y);
       // f = rotateZ(f, u_paramC.z);
        //sphereFold(f, u_paramC.w); // Sphere Inversion
        f *= c;
        v *= c;
        f += fc;

        orbit += length(start - f);
    }

    float z = length(f.xy) - fu;
    
    //fd +=  (f.y - u_paramD.w) * u_paramD.z;
    float d = 0.1 + max(z, abs(length(f.xy) * f.z) / sqrt(dot(f, f))) / abs(v);
    
    return float2(d, orbit);
}
*/


float udBox(float3 p, float3 b)
{
    return length(max(abs(p) - b, 0.0));
}

float sdBox(float3 p, float3 b)
{
    float3 d = abs(p) - b;
    return min(max(d.x, max(d.y, d.z)), 0.0) + length(max(d, 0.0));
}

float2 polycrust(float3 p)
{
    float3 dim = u_paramB.xyz;

    float3 t = u_paramA.xyz;
    
    float d = 1e10;

    float scale = u_paramC.x;

    float s = 1.0;
    float orbit = 0.0;
        
    float3 p0 = p;

    for (int i = 0; i < _LEVELS; i++)
    {
        p = rotate_vector(p - t / s, fromtwovectors(float3(0,1,0), u_paramD.xyz));

        d = min(d, sdBox(p.xyz / s, dim) * s);
        p = abs(p);

        // float circle =  fu/10.0 + 0.1 * sin(_Time.x + p.xyz);
        // d = min(d, length(p - t) -circle);

        s *= scale;
        orbit += length(p - p0);
        
    }

    return float2(d, orbit);
}

void boxFold(inout float3 z, inout float dz)
{
    z = clamp(z, -u_paramA.z, u_paramA.z) * 2.0 - z;
}

//https://www.shadertoy.com/view/4ds3zn
float2 appolonian( float3 p)
{
    float s = u_paramA.y;
	float scale = u_paramA.x;

	float4 orb = (float4)1000.0; 
	
	for( int i=0; i<_LEVELS ;i++ )
	{
		p = -1.0 + 2.0*frac(0.5*p+0.5);
        //p = rotate_vector(p, fromtwovectors(float3(0,1,0), u_paramB.xyz));

		float r2 = dot(p,p);
		
        orb = min( orb, float4(abs(p),r2) );
		
		float k = s/r2;
		p     *= k;
		scale *= k;
	}
	
	return float2(0.25*length(p)/scale, length(orb.y));
}


//https://www.shadertoy.com/view/4s3GW2
float2 Hoskins( float3 p )
{
	float scale = 1.;
    p= abs(p);
	float4 orb = (float4)1000.0; 

	for( int i=0; i < _LEVELS;i++ )
	{
		p = 2.0*clamp(p, -u_paramA.x, u_paramA.x) - p;

		float r2 = dot(p,p);
        //float r2 = dot(p,p+sin(p.z*.3)); //Alternate fractal

		float k = max(u_paramA.w/r2, u_paramA.y);
        orb = min( orb, float4(abs(p),r2) );

		p     *= k;
		scale *= k;
	}


	float l = length(p.xz);
	float rxy = l - u_paramA.z;
	float n = l * p.z;
	rxy = max(rxy, -(n) / u_paramA.y);
	
    
    return float2((rxy) / abs(scale), length(orb));
}


//----------------------------------------------------------------------------------------
float2 MBOX(float3 z)
{
    float3 offset = z;
    float dr = 1;

    float Scale = u_paramA.w;
    float iter = 0.0;

    float orbit = 0;
    float3 z_prime = z;

    for (int n = 0; n < _LEVELS; n++)
    {
        boxFold(z, dr); // Reflect
        sphereFold(z, dr); // Sphere Inversion
 		
        z = Scale * z + offset; // Scale & Translate
        dr = dr * abs(Scale) + 1.0;
        iter++;
        orbit += length(z_prime - z);
        z_prime = z;

        if (abs(dr) > 100000.)
            break;
    }
     
    float r = length(z);

    return float2(r / abs(dr), orbit);
}



//--------------------------------------------------------------------------------
// quaternion manipulation
//--------------------------------------------------------------------------------

float4 qSquare(float4 a)
{
    return float4(a.x * a.x - dot(a.yzw, a.yzw), 2.0 * a.x * (a.yzw));
}

float4 qCube(float4 a)
{
    return a * (4.0 * a.x * a.x - dot(a, a) * float4(3.0, 1.0, 1.0, 1.0));
}

//--------------------------------------------------------------------------------

float lengthSquared(float4 z) { return dot(z, z); }

float3 julia(float3 p)
{
    float4 c = _Quaternion;

    float4 z = float4(p, 0.2);

    float m2 = 0.0;
    float2 t = (float2)1e10;

    float dz2 = 1.0;
    for (int i = 0; i < _LEVELS; i++)
    {
        // |dz|² = |3z²|²
        dz2 *= 9.0 * lengthSquared(qSquare(z));

        // z = z^3 + c		
        z = qCube(z) + c;

        // stop under divergence		
        m2 = dot(z, z);
        if (m2 > 10000.0) break;

        // orbit trapping ( |z|² and z_x  )
        t = min(t, float2(m2, abs(z.x)));

    }

    // distance estimator: d(z) = 0.5·log|z|·|z|/|dz|   (see http://iquilezles.org/www/articles/distancefractals/distancefractals.htm)
    float d = 0.25 * log(m2) * sqrt(m2 / dz2);


    //d = max(d, udBox(p - u_paramC.xyz, u_paramD.xyz));

    return float3(d, t);
}
float qLength2(in float4 q) { return dot(q, q); }

float juliaTrap(in float3 p)
{
    float4 z = float4(p, 0.0);
    float dz2 = 1.0;
    float m2 = 0.0;
    float n = 0.0;
    float o = 1e10;

    for (int i = 0; i < _LEVELS; i++)
    {
        // z' = 3z² -> |z'|² = 9|z²|²
        dz2 *= 9.0 * qLength2(qSquare(z));

        // z = z³ + c		
        z = qCube(z) + u_paramB;

        // stop under divergence		
        m2 = qLength2(z);

        // orbit trapping : https://iquilezles.org/www/articles/orbittraps3d/orbittraps3d.htm
        o = min(o, length(z.xz - float2(0.45, 0.55)) - 0.1);

        // exit condition
        if (m2 > 256.0) break;
        n += 1.0;
    }

    // sdf(z) = log|z|·|z|/|dz| : https://iquilezles.org/www/articles/distancefractals/distancefractals.htm
    float d = 0.25 * log(m2) * sqrt(m2 / dz2);

    d = min(o, d);

    return float2(d, n);
}




// greetz 2 Mikael Hvidtfeldt Christensen
// http://blog.hvidtfeldts.net/index.php/2011/09/distance-estimated-3d-fractals-v-the-mandelbulb-different-de-approximations/
float2 mandelbulb(float3 pos) {

    float Power = u_paramA.x;
    float Bailout = u_paramA.y;

    float3 z = pos; //
    float dr = 1.0;
    float r = 0.0;
    float iter;

    for (int i = 0; i < _LEVELS; i++) {
        r = length(z);
        if (r > Bailout) break;

        // convert to polar coordinates
        float theta = acos(z.z / r) + u_paramA.z;
        float phi = atan2(z.y, z.x) + u_paramA.w;
        dr = pow(r, Power - 1.0) * Power * dr + 1.0;

        // scale and rotate the point
        float zr = pow(r, Power);
        theta = theta * Power;
        phi = phi * Power;

        // convert back to cartesian coordinates
        z = zr * float3(sin(theta) * cos(phi), sin(phi) * sin(theta), cos(theta));
        z += pos;
        iter++;
    }


    return float2( 0.5 * log(r) * r / dr, iter);
}


float2 rot2D (float2 q, float a)
{
  return q * cos (a) + q.yx * sin (a) * float2 (-1., 1.);
}

float2 IFSOLDERBOXY(float3 p) {
    float s = u_paramA.x ;

    float d = 1e5;
    float orbit = 0.0;
    float dp = d; 
    float3 offs = u_paramB.xyz;

    float amp = 1./s; // Analogous to layer amplitude.
    float level = -1;

    for(int i=0; i<_LEVELS; i++){
        float p0 = p;
        p.xy = rot2D(p.xy, u_paramC.z);
        p.yz = rot2D(p.yz, u_paramC.x);
        p.zx = rot2D(p.zx, u_paramC.y);

        p = abs(p);

        //mirrors
        p.xy += step(p.x, p.y)*(p.yx - p.xy);
        p.xz += step(p.x, p.z)*(p.zx - p.xz);
        p.yz += step(p.y, p.z)*(p.zy - p.yz);
        p=abs(p);

        // Stretching about an offset.
        p = p*s + offs*(1. - s);
        p -= step(p, offs*(1. - s)*.5)*offs*(1. - s);

        // p=abs(p);

        d = max(-d, max(max(p.x, p.y), p.z)*amp);

        // if ( dp != d){
        //     level = i;
        //     orbit ++;
        // }

        orbit += max(p.x, max(p.y, p.z));

        dp = d;
        amp /= s; 
    }

    return float2(d -  (u_paramA.z) * s , orbit);

}

float sdTorus( float3 p, float2 t )
{
  float2 q = float2(length(p.xz)-t.x,p.y);
  return length(q)-t.y;
}

float2 IFS_LEVELS(float3 p){
    float3 t = u_paramB.xyz;
    float3 offs =  u_paramD.xyz;
    float d = 1e10;
    float s = 1.0;
    float scale = u_paramA.x;

    p += t; 
    float level = 0;

    for (int i = 0; i < _LEVELS; i++){
  
        p = p-t/s;
        p = rotateX(p, u_paramC.x);
        p = rotateY(p, u_paramC.y);
        p = rotateZ(p, u_paramC.z);

        p.x = abs(p.x);
        p.y = abs(p.y);
        p.z = abs(p.z);

        float oldD = d;

        // box
       d = min	(d, udBox(p.xyz*s, offs)/s ) ;

        // torus
        // d = min	(d, sdTorus(p.xyz*s, offs.xy)/s ) ;
        
        // circle
        //d = min(d, length(p - t) - u_paramA.y);

        if ( d != oldD ) {
            level = i;
        }

        s *= scale * scale;
    }

    if ( level ==  round(u_paramA.z) ||  level ==  round(u_paramA.w)) {
        level = -1;
    }

    return float2(d, level);
}




float2 IFS(float3 p)
{
    float3 z = p;
    float3 Offset = u_paramB;
    float Scale = u_paramA.x;

    float r = 0;
    int n = 0;
    float dd = 99999;
    float level = -1;
    for(int i = 0; i < _LEVELS; i ++) {
        z = abs(z);
        z = rotateX(z, u_paramC.x);
        z = rotateY(z, u_paramC.y);
        z = rotateZ(z, u_paramC.z);
        z = z * Scale - Offset * (Scale - 1.0);
        n++; 
        
        float oldD = dd;
        dd = min(dd, sdBox(z, u_paramD.xyz));
        if ( dd != oldD ) {
            level = n;
        }

    }

    /* if (sDist < fDist)
     {
         return float2(sDist, 5);
     }
     else*/
 
    if (round(u_paramA.y) == level) {
        level = -1;
    } 

    return float2(dd, level);                                                                                                                                           
}

float opSmoothIntersection(float d1, float d2, float k) {
    float h = clamp(0.5 - 0.5 * (d2 - d1) / k, 0.0, 1.0);
    return lerp(d2, d1, h) + k * h * (1.0 - h);
}

// float2 DE(float3 p)
// {

// //     float2 f = appolonian(p);
// //    float d =  sdBox(p- u_paramC.xyz, float3(u_paramD.x, u_paramD.y, u_paramD.z));
//     // float d =  length(p - u_paramC.xyz ) - u_paramD.x;
    
// //     if ( min(f.x, d) == f.x) {
// //         return f;
// //     } else {
// //         // -1 indicates emissive
// //         return float2(d, -1);
// //     }



//     float2 d = pseudo_knightyan(p);


//     // return float2(d, f.y);
//     //float2 d = juliaTrap(p);
//    float d2 = length(p) - u_paramD.x; // + d.y * u_paramD.y;

//    return float2(max(d.x, d2), d.y);
//   //  return min(tglad_variant(p),);
//     // return tglad_variant(p);
// }

float2 DE(float3 p)
{

    float2 d = float2(0,0);
    
    if ( _FractalType < 1) {
        d = pseudo_knightyan(p);

    } else if (_FractalType < 2) {
        d = tglad_variant(p);

    } else if (_FractalType < 3) {

        d = tglad(p);

    } else if (_FractalType < 4) {

        d = MBOX(p);

    } else if (_FractalType < 5) {

        d = mandelbulb(p);

    } else if (_FractalType < 6) {

        d = hartverdrahtet(p);

    } else if (_FractalType < 7) {

        d = appolonian(p);

    } else if (_FractalType < 8) {

        d = juliaTrap(p);

    }
//    float d2 =  sdBox(p- u_paramC.xyz, float3(u_paramD.x, u_paramD.y, u_paramD.z));
    // float d =  length(p - u_paramC.xyz ) - u_paramD.x;
    
//     if ( min(f.x, d) == f.x) {
//         return f;
//     } else {
//         // -1 indicates emissive
//         return float2(d, -1);
//     }


    // float2 d = tglad_variant(p);


    // return float2(d, f.y);
//    float d2 = length(p) - u_paramD.x; // + d.y * u_paramD.y;

//    return float2(max(d.x, d2), d.y);
  //  return min(tglad_variant(p),);
//    return float2(min(d.x, d2), d.y);
    return d;

}