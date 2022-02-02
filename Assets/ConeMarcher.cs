//using FFmpegOut;
using System.Collections;
using System.Collections.Generic;
using UnityEditor;
using UnityEngine;
using System.IO;
using System;
using UnityEngine.Playables;

[ExecuteInEditMode]
public class ConeMarcher : MonoBehaviour
{
    public ComputeShader ConeMarchingKernel;
    public ComputeShader ShadingKernel;
    public ComputeShader ReconstructionKernel;
    public ComputeShader RefineKernel;
 
    public Texture _matcapTexture;
    private Camera _camera;

    public Light _light;

    RenderTexture[] _hierarchy;
    RenderTexture _reconstruction;
    RenderTexture _currentShadedFrame;

    public float Threshold = 0.1f;

    private int _currentSample = 0;

    // level 0 is fulll res
    private int NUM_LEVELS = 8;

    public Vector4 _ParamA;
    public Vector4 _ParamB;
    public Vector4 _ParamC;
    public Vector4 _ParamD;

    private void Awake()
    {
        _camera = GetComponent<Camera>();
    }

    public void SetDirty()
    {

    }

    private void OnEnable()
    {
    }
     
    private void Update()
    {
     
    }

    Matrix4x4 _lastCamToWorld;

    private void SetShaderParameters()
    {
        // var camToWorld = _camera.transform.localToWorldMatrix;
        var camToWorld = _camera.cameraToWorldMatrix;
        var invProjection = _camera.projectionMatrix.inverse;

        var pixelwidth = 2 * Mathf.Atan(0.5f * Mathf.Deg2Rad * _camera.fieldOfView /Screen.width) ; 

        ConeMarchingKernel.SetFloat("_PixelWidth", pixelwidth);
        ConeMarchingKernel.SetMatrix("_CameraToWorld", camToWorld);
        ConeMarchingKernel.SetMatrix("_InverseProjection", invProjection);
        ConeMarchingKernel.SetVector("_ClipPlanes", new Vector2(_camera.nearClipPlane, _camera.farClipPlane));
        ConeMarchingKernel.SetVector("_ParamA", _ParamA);
        ConeMarchingKernel.SetVector("_ParamB", _ParamB);
        ConeMarchingKernel.SetVector("_ParamC", _ParamC);
        ConeMarchingKernel.SetVector("_ParamD", _ParamD);

        // RefineKernel.SetTexture(0, "_Previous", _hierarchy[0]);
        RefineKernel.SetTexture(0, "_Result", _hierarchy[0]);
        RefineKernel.SetInt("_Level", 0);
        RefineKernel.SetVector("_Resolution", new Vector2(_hierarchy[0].width, _hierarchy[0].height));
        RefineKernel.SetVector("_LightDir", _light.transform.forward);
        RefineKernel.SetFloat("_PixelWidth", pixelwidth);
        RefineKernel.SetMatrix("_CameraToWorld", camToWorld);
        RefineKernel.SetMatrix("_InverseProjection", invProjection);
        RefineKernel.SetVector("_ParamA", _ParamA);
        RefineKernel.SetVector("_ParamB", _ParamB);
        RefineKernel.SetVector("_ParamC", _ParamC);
        RefineKernel.SetVector("_ParamD", _ParamD);

        ShadingKernel.SetTexture(0, "_Matcap", _matcapTexture);
        ShadingKernel.SetTexture(0, "_GBuffer", _hierarchy[0]);
        ShadingKernel.SetTexture(0, "_Result", _currentShadedFrame);
        ShadingKernel.SetInt("_Level", 0);
        ShadingKernel.SetVector("_Resolution", new Vector2(_currentShadedFrame.width, _currentShadedFrame.height));
        ShadingKernel.SetVector("_LightDir", _light.transform.forward);
        ShadingKernel.SetFloat("_PixelWidth", pixelwidth);
        ShadingKernel.SetMatrix("_CameraToWorld", camToWorld);
        ShadingKernel.SetMatrix("_InverseProjection", invProjection);
        ShadingKernel.SetVector("_ParamA", _ParamA);
        ShadingKernel.SetVector("_ParamB", _ParamB);
        ShadingKernel.SetVector("_ParamC", _ParamC);
        ShadingKernel.SetVector("_ParamD", _ParamD);

        ReconstructionKernel.SetTexture(0, "_Current", _currentShadedFrame);
        ReconstructionKernel.SetTexture(0, "_Reconstruction", _reconstruction);
        ReconstructionKernel.SetMatrix("_LastFrameCameraToWorld", _lastCamToWorld);
        ReconstructionKernel.SetMatrix("_CameraToWorld", camToWorld);
        ReconstructionKernel.SetMatrix("_InverseProjection", invProjection);

    }

    private RenderTexture createTexture(int w, int h) {
        RenderTexture r = new RenderTexture(Screen.width, Screen.height, 0, RenderTextureFormat.ARGBFloat, RenderTextureReadWrite.Linear);
        r.enableRandomWrite = true;
        r.Create();
        return r;
    }

    private void CreateTextures()
    {
        if (_reconstruction == null || _reconstruction.width != Screen.width || _reconstruction.height != Screen.height)
        {

            _reconstruction?.Release();
            _currentShadedFrame?.Release();

            if (_hierarchy != null){
                for( int i = 0 ; i < _hierarchy.Length; i++) {
                    _hierarchy[i]?.Release();
                }
            }
 
            _hierarchy = new RenderTexture[NUM_LEVELS];
            for(int i = 0; i< NUM_LEVELS; i ++) {
                float multiplier = 1f / Mathf.Pow(2, i);
                int w = Mathf.CeilToInt(Screen.width * multiplier);
                int h = Mathf.CeilToInt(Screen.height * multiplier);
                _hierarchy[i] = createTexture(w,h);
            }

            int sw = Screen.width;
            int sh = Screen.height;

            _currentShadedFrame = createTexture(sw,sh);
            _reconstruction = createTexture(sw,sh);
            
        }
    }

    private void Render(RenderTexture destination)
    {
        CreateTextures();
        SetShaderParameters();

        for(int i = NUM_LEVELS-1; i >= 0; i--) {
           
            if ( i < NUM_LEVELS-1) {
                ConeMarchingKernel.SetTexture(0, "_Previous", _hierarchy[i+1] );
            }

            ConeMarchingKernel.SetInt("_Level", i);
            ConeMarchingKernel.SetTexture(0, "_Result", _hierarchy[i]);
            ConeMarchingKernel.SetVector("_Resolution", new Vector2(_hierarchy[i].width, _hierarchy[i].height));

            int threadGroupsX = Mathf.CeilToInt(_hierarchy[i].width / 8.0f);
            int threadGroupsY = Mathf.CeilToInt(_hierarchy[i].height / 8.0f);
            ConeMarchingKernel.Dispatch(0, threadGroupsX, threadGroupsY, 1);
        }

        int tx = Mathf.CeilToInt(_hierarchy[0].width / 8.0f);
        int ty = Mathf.CeilToInt(_hierarchy[0].height / 8.0f);

        // RefineKernel.Dispatch(0, tx, ty, 1);
        ShadingKernel.Dispatch(0, tx, ty, 1);
        ReconstructionKernel.Dispatch(0, tx, ty, 1);

        Graphics.Blit(_reconstruction, destination);
    }

    private void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        Render(destination);
        _currentSample++;
        _lastCamToWorld = _camera.worldToCameraMatrix;
    }

    private void OnValidate()
    {
        _currentSample = 0;
    }
}
