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
    
    private ComputeBuffer _sphereBuffer;
    
    public static readonly int SphereBufferPropertyID = Shader.PropertyToID("SphereBuffer");
    public static readonly int SphereCountPropertyID = Shader.PropertyToID("SphereCount");
    public static readonly int CameraFocalDistancePropertyID = Shader.PropertyToID("CameraFocalDistance");
    public static readonly int CameraPlaneWidthPropertyID = Shader.PropertyToID("CameraPlaneWidth");
    public static readonly int CameraPlaneHeightPropertyID = Shader.PropertyToID("CameraPlaneHeight");
    public static readonly int SamplesPerPixelPropertyID = Shader.PropertyToID("SamplesPerPixel");
    public static readonly int PixelDeltaUPropertyID = Shader.PropertyToID("PixelDeltaU");
    public static readonly int PixelDeltaVPropertyID = Shader.PropertyToID("PixelDeltaV");
    public static readonly int Pixel00LocPropertyID = Shader.PropertyToID("Pixel00Loc");
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
        float focalDistance = 1.0f;
        
        // Calculate the dimensions of the camera's near plane
        float planeHeight = 2.0f * focalDistance * Mathf.Tan(camera.fieldOfView * 0.5f * Mathf.Deg2Rad);
        float planeWidth = planeHeight * camera.aspect;
        
        // Calculate pixel spacing vectors for antialiasing
        Vector3 pixelDeltaU = new Vector3(planeWidth / Screen.width, 0, 0);
        Vector3 pixelDeltaV = new Vector3(0, -planeHeight / Screen.height, 0); // Negative because screen Y is flipped
        
        // Calculate the location of pixel (0,0) in world space relative to camera
        Vector3 viewportUpperLeft = new Vector3(-planeWidth * 0.5f, planeHeight * 0.5f, focalDistance);
        Vector3 pixel00Loc = viewportUpperLeft + 0.5f * (pixelDeltaU + pixelDeltaV);
        
        // Populate camera data
        material.SetFloat(CameraFocalDistancePropertyID, focalDistance);
        material.SetFloat(CameraPlaneWidthPropertyID, planeWidth);
        material.SetFloat(CameraPlaneHeightPropertyID, planeHeight);
        material.SetVector(PixelDeltaUPropertyID, pixelDeltaU);
        material.SetVector(PixelDeltaVPropertyID, pixelDeltaV);
        material.SetVector(Pixel00LocPropertyID, pixel00Loc);
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