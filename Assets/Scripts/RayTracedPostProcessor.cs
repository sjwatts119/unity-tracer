using System;
using Core;
using Geometry;
using UnityEngine;
using UnityEngine.Internal;


[ExecuteAlways, ImageEffectAllowedInSceneView]
public class RayTracedPostProcessor : MonoBehaviour
{
    [Header("Anti-Aliasing")]
    [Range(1, 250)]
    public int samplesPerPixel;
    
    [Header("Ray Tracing")]
    [Range(1, 250)]
    public int rayMaxDepth;
    
    [Header("Camera")]
    [Range(0.1f, 50f)]
    public float cameraFocalDistance = 1.0f;
    
    [Header("Camera")]
    [Range(0.0f, 100f)]
    public float cameraDefocusAngle = 10f;
    
    private ComputeBuffer _sphereBuffer;
    
    public static readonly int SphereBufferPropertyID = Shader.PropertyToID("SphereBuffer");
    public static readonly int SphereCountPropertyID = Shader.PropertyToID("SphereCount");
    public static readonly int CameraFocalDistancePropertyID = Shader.PropertyToID("CameraFocalDistance");
    public static readonly int CameraPlaneWidthPropertyID = Shader.PropertyToID("CameraPlaneWidth");
    public static readonly int CameraPlaneHeightPropertyID = Shader.PropertyToID("CameraPlaneHeight");
    public static readonly int SamplesPerPixelPropertyID = Shader.PropertyToID("SamplesPerPixel");
    public static readonly int RayMaxDepthPropertyID = Shader.PropertyToID("RayMaxDepth");

    // Called after all rendering is complete to render image
    public void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        var material = RayShader.Material;
        
        PopulateAntialiasingData(material);
        PopulateRaytracingData(material);
        PopulateSphereData(material);
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
        material.SetFloat("CameraDefocusAngle", cameraDefocusAngle);
        material.SetFloat(CameraPlaneWidthPropertyID, planeWidth);
        material.SetFloat(CameraPlaneHeightPropertyID, planeHeight);
        material.SetMatrix("CameraLocalToWorld", camera.transform.localToWorldMatrix);
        
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

    // Release the compute buffer when disabled or destroyed
    void OnDisable() => _sphereBuffer?.Release();
    void OnDestroy() => _sphereBuffer?.Release();
}