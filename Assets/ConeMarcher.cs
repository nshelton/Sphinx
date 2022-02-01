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
 
    private Camera _camera;

    RenderTexture[] _hierarchy;
    RenderTexture _reconstruction;

    public float Threshold = 0.1f;

    private int _currentSample = 0;

    // level 0 is fulll res
    private int NUM_LEVELS = 8;

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

    private void SetShaderParameters()
    {
        ConeMarchingKernel.SetFloat("_Threshold", Threshold);
        ConeMarchingKernel.SetMatrix("_CameraToWorld", _camera.cameraToWorldMatrix);

        ShadingKernel.SetTexture(0, "_GBuffer", _hierarchy[0]);
        ShadingKernel.SetTexture(0, "_Result", _reconstruction);
        ShadingKernel.SetInt("_Level", 0);
        ShadingKernel.SetVector("_Resolution", new Vector2(_reconstruction.width, _reconstruction.height));
        ShadingKernel.SetFloat("_Threshold", Threshold);
        ShadingKernel.SetMatrix("_CameraToWorld", _camera.cameraToWorldMatrix);
    }

    private void CreateTextures()
    {
        if (_reconstruction == null || _reconstruction.width != Screen.width || _reconstruction.height != Screen.height)
        {
            // Release render texture if we already have one
            if (_reconstruction != null)
            {
                _reconstruction.Release();
            }

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
                _hierarchy[i] = new RenderTexture(w, h, 0, RenderTextureFormat.ARGBFloat, RenderTextureReadWrite.Linear);
                _hierarchy[i].enableRandomWrite = true;
                _hierarchy[i].Create();
            }

            _reconstruction = new RenderTexture(Screen.width, Screen.height, 0, RenderTextureFormat.ARGBFloat, RenderTextureReadWrite.Linear);
            _reconstruction.enableRandomWrite = true;
            _reconstruction.Create();
        }
    }

    public void ClearRendertexture(RenderTexture renderTexture)
    {
        RenderTexture rt = RenderTexture.active;
        RenderTexture.active = renderTexture;
        GL.Clear(true, true, Color.clear);
        RenderTexture.active = rt;
    }

    private void Render(RenderTexture destination)
    {
        // Make sure we have a current render target
        CreateTextures();
        SetShaderParameters();

        // Set the target and dispatch the compute shader
 
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

        ShadingKernel.Dispatch(0, tx, ty, 1);

        Graphics.Blit(_reconstruction, destination);
    }

    Matrix4x4 m_worldToLastFrame;

    private void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        Render(destination);
        _currentSample++;
        m_worldToLastFrame = _camera.worldToCameraMatrix;
    }

    private void OnValidate()
    {
        _currentSample = 0;
    }
}
