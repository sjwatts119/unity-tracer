using System;
using Core;
using Geometry;
using UnityEngine;

[ExecuteAlways, ImageEffectAllowedInSceneView]
public class RayTracedPostProcessor : MonoBehaviour
{
    private ComputeBuffer _sphereBuffer;
    private static readonly int SphereBufferPropertyID = Shader.PropertyToID("SphereBuffer");

    // Called after all rendering is complete to render image
    public void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        var material = RayShader.Material;
        
        PopulateSphereBuffer(material);
        
        Graphics.Blit(source, destination, material);
    }

    // Populate the gpu compute buffer with sphere data from the scene
    void PopulateSphereBuffer(Material material)
    {
        var spheres = FindObjectsByType<Sphere>(FindObjectsSortMode.None);
        
        _sphereBuffer?.Release();
        _sphereBuffer = null;
        
        if (spheres.Length > 0)
        {
            var sphereData = new ShaderStructs.Sphere[spheres.Length];
            for (int i = 0; i < spheres.Length; i++)
            {
                sphereData[i] = spheres[i].ToShaderData();
            }

            _sphereBuffer = new ComputeBuffer(sphereData.Length, System.Runtime.InteropServices.Marshal.SizeOf<ShaderStructs.Sphere>());
            _sphereBuffer.SetData(sphereData);
            material.SetBuffer(SphereBufferPropertyID, _sphereBuffer);
        }
    }

    // Release the compute buffer when disabled or destroyed
    void OnDisable() => _sphereBuffer?.Release();
    void OnDestroy() => _sphereBuffer?.Release();
}