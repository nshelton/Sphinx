    float4 _ParamA;
float4 _ParamB;
float4 _ParamC;
float4 _ParamD;

float tglad_formula(float3 z0)
{
    float sp = length(z0+float3(0.0, 2.0, 0.0));
    z0 = frac(z0 * 2.0 ) * 2.0;
    // z0 = frac(z0 * 2.0 ) * 2.0;

    float mr=0.25, mxr=1.0;
    float4 scale=float4(-3.12,-3.12,-3.12,3.12);
    float4 p0=float4(0.0,1.59,-1.0,0.0);
    float4 z = float4(z0,1.0);

    for (int n = 0; n < 3; n++) {
        z.xyz=clamp(z.xyz, -0.94, 0.94)*2.0-z.xyz;
        z*=scale/clamp(dot(z.xyz,z.xyz),mr,mxr);
        z+=p0;
    }

    float dS=(length(max(abs(z.xyz)-float3(1.2,49.0,1.4),0.0))-0.06)/z.w;
    dS = max(sp, dS);
    return dS;
}

void boxFold(inout float3 z, inout float dz)
{
    z = clamp(z, -_ParamA.z, _ParamA.z) * 2.0 - z;
}

void sphereFold(inout float3 z, inout float dz)
{
    float fixedRadius2 = _ParamA.x;
    float minRadius2 = _ParamA.y;

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


//https://www.shadertoy.com/view/4ds3zn
float appolonian( float3 p)
{
    float s = _ParamA.y;
	float scale = _ParamA.x;

	float4 orb = (float4)1000.0; 
	
	for( int i=0; i<6 ;i++ )
	{
		p = -1.0 + 2.0*frac(0.5*p+0.5);
        //p = rotate_vector(p, fromtwovectors(float3(0,1,0), _ParamB.xyz));

		float r2 = dot(p,p);
		
        orb = min( orb, float4(abs(p),r2) );
		
		float k = s/r2;
		p     *= k;
		scale *= k;
	}
	// return float2(0.25*length(p)/scale, length(orb.y));
	
	return 0.25*length(p)/scale;
}

float MBOX(float3 z)
{
    float3 offset = z;
    float dr = 1;

    float Scale = _ParamA.w;
    float iter = 0.0;

    float orbit = 0;
    float3 z_prime = z;

    for (int n = 0; n < 12; n++)
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

    // return float2(r / abs(dr), orbit);
    return r / abs(dr);
}


float hartverdrahtet(float3 f, int i)
{
    float3 cs=float3(.808,.808,1.167);
    float fs=1.;
    float3 fc=0;
    float fu=10.;
    float fd=.763;
    
    // scene selection
    {
        if(i==0) cs.y=.58;
        if(i==1) cs.xy=.5;
        if(i==2) cs.xy=.5;
        if(i==3) fu=1.01,cs.x=.9;
        if(i==4) fu=1.01,cs.x=.9;
        if(i==6) cs=float3(.5,.5,1.04);
        if(i==5) fu=.9;
        if(i==7) fd=.7,fs=1.34,cs.xy=.5;
        if(i==8) fc.z=-.38;
    }
    
    //cs += sin(time)*0.2;

    float v=1.;
    for(int i=0; i<12; i++){
        f=2.*clamp(f,-cs,cs)-f;
        float c=max(fs/dot(f,f),1.);
        f*=c;
        v*=c;
        f+=fc;
    }
    float z=length(f.xy)-fu;
    return fd*max(z,abs(length(f.xy)*f.z)/sqrt(dot(f,f)))/abs(v);
}

float map(float3 p ){
    return MBOX(p);
    // return length(p) - 0.5;
}