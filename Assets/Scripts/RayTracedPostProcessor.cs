using System;
using System.Collections.Generic;
using Core;
using Geometry;
using NUnit.Framework;
using Unity.VisualScripting;
using UnityEngine;
using Plane = Geometry.Plane;


[ExecuteAlways, ImageEffectAllowedInSceneView]
public class RayTracedPostProcessor : MonoBehaviour
{
    [Header("Anti-Aliasing")]
    [UnityEngine.Range(1, 250)]
    public int samplesPerPixel;
    
    [Header("Ray Tracing")]
    [UnityEngine.Range(1, 250)]
    public int rayMaxDepth;
    
    [Header("Camera")]
    [UnityEngine.Range(0.1f, 50f)]
    public float cameraFocalDistance = 1.0f;
    
    [Header("Camera")]
    [UnityEngine.Range(0.0f, 100f)]
    public float cameraDefocusAngle = 10f;
    
    private ComputeBuffer _sphereBuffer;
    private ComputeBuffer _quadBuffer;
    
    public static readonly int SphereBufferPropertyID = Shader.PropertyToID("SphereBuffer");
    public static readonly int SphereCountPropertyID = Shader.PropertyToID("SphereCount");
    public static readonly int QuadBufferPropertyID = Shader.PropertyToID("QuadBuffer");
    public static readonly int QuadCountPropertyID = Shader.PropertyToID("QuadCount");
    public static readonly int CameraFocalDistancePropertyID = Shader.PropertyToID("CameraFocalDistance");
    public static readonly int CameraPlaneWidthPropertyID = Shader.PropertyToID("CameraPlaneWidth");
    public static readonly int CameraPlaneHeightPropertyID = Shader.PropertyToID("CameraPlaneHeight");
    public static readonly int CameraDefocusAnglePropertyID = Shader.PropertyToID("CameraDefocusAngle");
    public static readonly int CameraLocalToWorldPropertyID = Shader.PropertyToID("CameraLocalToWorld");
    public static readonly int SamplesPerPixelPropertyID = Shader.PropertyToID("SamplesPerPixel");
    public static readonly int RayMaxDepthPropertyID = Shader.PropertyToID("RayMaxDepth");

    // Called after all rendering is complete to render image
    public void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        var material = RayShader.Material;
        
        PopulateAntialiasingData(material);
        PopulateRaytracingData(material);
        PopulateSphereData(material);
        PopulateQuadData(material);
        PopulateCameraData(material, GetComponent<Camera>());
        
        Graphics.Blit(source, destination, material);
    }
    
    void PopulateCameraData(Material material, Camera camera)
    {
        // Calculate the dimensions of the camera's near plane
        float planeHeight = 2.0f * cameraFocalDistance * Mathf.Tan(camera.fieldOfView * 0.5f * Mathf.Deg2Rad);
        float planeWidth = planeHeight * camera.aspect;
        
        // Populate camera data
        material.SetFloat(CameraFocalDistancePropertyID, cameraFocalDistance);
        material.SetFloat(CameraDefocusAnglePropertyID, cameraDefocusAngle);
        material.SetFloat(CameraPlaneWidthPropertyID, planeWidth);
        material.SetFloat(CameraPlaneHeightPropertyID, planeHeight);
        material.SetMatrix(CameraLocalToWorldPropertyID, camera.transform.localToWorldMatrix);
    }
    
    void PopulateAntialiasingData(Material material)
    {
        material.SetInt(SamplesPerPixelPropertyID, samplesPerPixel);
    }
    
    void PopulateRaytracingData(Material material)
    {
        material.SetInt(RayMaxDepthPropertyID, rayMaxDepth);
    }

    // Populate the gpu compute buffer with sphere data from the scene
    void PopulateSphereData(Material material)
    {
        var spheres = FindObjectsByType<Sphere>(FindObjectsSortMode.None);
        
        _sphereBuffer?.Release();
        _sphereBuffer = null;
        
        if (spheres.Length == 0)
        {
            material.SetInt(SphereCountPropertyID, 0);
            return;
        }

        var sphereData = new ShaderStructs.Sphere[spheres.Length];
        for (int i = 0; i < spheres.Length; i++)
        {
            sphereData[i] = spheres[i].ToShaderData();
        }

        // Populate the compute buffer with sphere data
        _sphereBuffer = new ComputeBuffer(sphereData.Length, System.Runtime.InteropServices.Marshal.SizeOf<ShaderStructs.Sphere>());
        _sphereBuffer.SetData(sphereData);
        
        // Bind the buffer and count to the shader
        material.SetBuffer(SphereBufferPropertyID, _sphereBuffer);
        material.SetInt(SphereCountPropertyID, sphereData.Length);
    }

// Replace the PopulateQuadData method in your RayTracedPostProcessor class
    void PopulateQuadData(Material material)
    {
        // Find both quads and planes
        var quads = FindObjectsByType<Quad>(FindObjectsSortMode.None);
        var planes = FindObjectsByType<Plane>(FindObjectsSortMode.None);
    
        _quadBuffer?.Release();
        _quadBuffer = null;
    
        int totalCount = quads.Length + planes.Length;
    
        if (totalCount == 0)
        {
            material.SetInt(QuadCountPropertyID, 0);
            return;
        }
    
        var quadData = new ShaderStructs.Quad[totalCount];
        int index = 0;
    
        for (int i = 0; i < quads.Length; i++)
        {
            quadData[index++] = quads[i].ToShaderData();
        }
    
        // Add all planes to the buffer
        for (int i = 0; i < planes.Length; i++)
        {
            quadData[index++] = planes[i].ToShaderData();
        }
    
        // Populate the compute buffer with combined quad/plane data
        _quadBuffer = new ComputeBuffer(quadData.Length, System.Runtime.InteropServices.Marshal.SizeOf<ShaderStructs.Quad>());
        _quadBuffer.SetData(quadData);
    
        // Bind the buffer and count to the shader
        material.SetBuffer(QuadBufferPropertyID, _quadBuffer);
        material.SetInt(QuadCountPropertyID, quadData.Length);
    }

    // Release the compute buffer when disabled or destroyed
    void OnDisable()
    {
        _sphereBuffer?.Release();
        _quadBuffer?.Release();
    }

    void OnDestroy()
    {
        _sphereBuffer?.Release();
        _quadBuffer?.Release();
    }
}