using System;
using Core;
using Geometry;
using UnityEngine;


[ExecuteAlways, ImageEffectAllowedInSceneView]
public class RayTracedPostProcessor : MonoBehaviour
{
    private ComputeBuffer _sphereBuffer;
    public static readonly int SphereBufferPropertyID = Shader.PropertyToID("SphereBuffer");
    public static readonly int SphereCountPropertyID = Shader.PropertyToID("SphereCount");
    public static readonly int CameraFocalDistancePropertyID = Shader.PropertyToID("CameraFocalDistance");
    public static readonly int CameraPlaneWidthPropertyID = Shader.PropertyToID("CameraPlaneWidth");
    public static readonly int CameraPlaneHeightPropertyID = Shader.PropertyToID("CameraPlaneHeight");
    

    // Called after all rendering is complete to render image
    public void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        var material = RayShader.Material;
        
        PopulateSphereBuffer(material);
        PopulateCameraDataBuffer(material, GetComponent<Camera>());
        
        Graphics.Blit(source, destination, material);
    }
    
    void PopulateCameraDataBuffer(Material material, Camera camera)
    {
        float focalDistance = 1.0f;
        
        // Calculate the dimensions of the camera's near plane
        float planeHeight = 2.0f * focalDistance * Mathf.Tan(camera.fieldOfView * 0.5f * Mathf.Deg2Rad);
        float planeWidth = planeHeight * camera.aspect;
        
        // Populate camera data
        material.SetFloat(CameraFocalDistancePropertyID, focalDistance);
        material.SetFloat(CameraPlaneWidthPropertyID, planeWidth);
        material.SetFloat(CameraPlaneHeightPropertyID, planeHeight);
    }

    // Populate the gpu compute buffer with sphere data from the scene
    void PopulateSphereBuffer(Material material)
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