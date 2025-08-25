using System;
using System.Linq;
using Core;
using Geometry;
using UnityEngine;
using Plane = Geometry.Plane;

[ExecuteAlways, ImageEffectAllowedInSceneView]
public class RayTracedPostProcessor : MonoBehaviour
{
    [Header("Info")] 
    public int renderedFrames;
    
    [Header("Anti-Aliasing")]
    [UnityEngine.Range(1, 1000)]
    public int samplesPerPixel;
    
    [Header("Ray Tracing")]
    [UnityEngine.Range(1, 250)]
    public int rayMaxDepth;
    
    [Header("Rendering")]
    public bool accumulateFrames = true;
    
    [Header("Camera")]
    [UnityEngine.Range(0.1f, 50f)]
    public float cameraFocalDistance = 1.0f;
    
    [UnityEngine.Range(0.0f, 100f)]
    public float cameraDefocusAngle = 10f;

    [Header("Accumulation")]
    [SerializeField] private Shader accumulationShader;
    
    private ComputeBuffer _sphereBuffer;
    private ComputeBuffer _quadBuffer;
    private RenderTexture _accumulatedTexture;
    private Material _accumulationMaterial;
    
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
    public static readonly int FrameNumberPropertyID = Shader.PropertyToID("FrameNumber");
    public static readonly int PreviousTexturePropertyID = Shader.PropertyToID("PreviousTexture");

    private void Start()
    {
        renderedFrames = 0;
        InitialiseAccumulation();
    }

    // Called after all rendering is complete to render image
    public void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        var rayMaterial = RayShader.Material;
        
        // Populate shader data from scene and settings
        PopulateAntialiasingData(rayMaterial);
        PopulateRaytracingData(rayMaterial);
        PopulateSphereData(rayMaterial);
        PopulateQuadData(rayMaterial);
        PopulateCameraData(rayMaterial, GetComponent<Camera>());
        PopulateRenderingData(rayMaterial);

        if (!accumulateFrames)
        {
            renderedFrames = 0;
            Graphics.Blit(source, destination, rayMaterial);
            return;
        }

        InitialiseAccumulation();
        CreateAccumulationTexture(source.width, source.height);

        // Create temporary copy of the previous accumulated frame
        RenderTexture prevFrameCopy = RenderTexture.GetTemporary(source.width, source.height, 0, RenderTextureFormat.ARGBFloat);
        Graphics.Blit(_accumulatedTexture, prevFrameCopy);

        // Render the current frame with ray tracing
        RenderTexture currentFrame = RenderTexture.GetTemporary(source.width, source.height, 0, RenderTextureFormat.ARGBFloat);
        Graphics.Blit(source, currentFrame, rayMaterial);

        // Blend current frame with accumulated previous frames
        if (_accumulationMaterial != null)
        {
            _accumulationMaterial.SetInt(FrameNumberPropertyID, renderedFrames);
            _accumulationMaterial.SetTexture(PreviousTexturePropertyID, prevFrameCopy);
            Graphics.Blit(currentFrame, _accumulatedTexture, _accumulationMaterial);
        }
        else
        {
            // just overwrite if no accumulation material is set
            Graphics.Blit(currentFrame, _accumulatedTexture);
        }

        // Present the accumulated result to screen
        Graphics.Blit(_accumulatedTexture, destination);

        // Clean up temporary textures
        RenderTexture.ReleaseTemporary(currentFrame);
        RenderTexture.ReleaseTemporary(prevFrameCopy);

        renderedFrames++;
    }

    void OnDisable()
    {
        ReleaseBuffers();
        ReleaseTextures();
    }

    void OnDestroy()
    {
        ReleaseBuffers();
        ReleaseTextures();
    }
    
    // Reset accumulation when parameters change
    void OnValidate()
    {
        if (Application.isPlaying)
        {
            renderedFrames = 0;
            ReleaseTextures();
        }
    }

    private void InitialiseAccumulation()
    {
        // Create accumulation material if needed
        if (accumulationShader != null && _accumulationMaterial == null)
        {
            _accumulationMaterial = new Material(accumulationShader);
        }
    }

    private void CreateAccumulationTexture(int width, int height)
    {
        // Create or recreate the accumulation texture if size has changed
        if (_accumulatedTexture == null || _accumulatedTexture.width != width || _accumulatedTexture.height != height)
        {
            if (_accumulatedTexture != null)
            {
                _accumulatedTexture.Release();
            }

            _accumulatedTexture = new RenderTexture(width, height, 0, RenderTextureFormat.ARGBFloat);
            _accumulatedTexture.enableRandomWrite = true;
            _accumulatedTexture.Create();
        }
    }

    /*
     * Passing data to shaders
     */
    
    void PopulateCameraData(Material material, Camera camera)
    {
        // Calculate the dimensions of the camera's near plane
        float planeHeight = 2.0f * cameraFocalDistance * Mathf.Tan(camera.fieldOfView * 0.5f * Mathf.Deg2Rad);
        float planeWidth = planeHeight * camera.aspect;
        
        // Send camera parameters to shader
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
    
    void PopulateRenderingData(Material material)
    {
        material.SetInt(FrameNumberPropertyID, accumulateFrames ? renderedFrames : 0);
    }

    void PopulateSphereData(Material material)
    {
        var spheres = FindObjectsByType<Sphere>(FindObjectsSortMode.None);
        
        // Clean up existing buffer
        _sphereBuffer?.Release();
        _sphereBuffer = null;
        
        // Calculate total count
        int totalCount = 0;
        foreach (var sphere in spheres) totalCount += sphere.ToShaderData().Length;

        // Early exit if no spheres in scene
        if (totalCount == 0)
        {
            material.SetInt(SphereCountPropertyID, 0);
            return;
        }

        // Collect all sphere data
        var allSphereData = new ShaderStructs.Sphere[totalCount];
        int index = 0;
        
        foreach (var sphere in spheres)
        {
            var data = sphere.ToShaderData();
            for (int i = 0; i < data.Length; i++)
                allSphereData[index++] = data[i];
        }

        // Create and populate compute buffer
        _sphereBuffer = new ComputeBuffer(allSphereData.Length, System.Runtime.InteropServices.Marshal.SizeOf<ShaderStructs.Sphere>());
        _sphereBuffer.SetData(allSphereData);
        
        // Bind buffer to shader
        material.SetBuffer(SphereBufferPropertyID, _sphereBuffer);
        material.SetInt(SphereCountPropertyID, allSphereData.Length);
    }

    void PopulateQuadData(Material material)
    {
        var quads = FindObjectsByType<Quad>(FindObjectsSortMode.None);
        var planes = FindObjectsByType<Plane>(FindObjectsSortMode.None);
        var cuboids = FindObjectsByType<Cuboid>(FindObjectsSortMode.None);

        _quadBuffer?.Release();
        _quadBuffer = null;

        // Calculate total count
        int totalCount = 0;
        foreach (var quad in quads) totalCount += quad.ToShaderData().Length;
        foreach (var plane in planes) totalCount += plane.ToShaderData().Length;
        foreach (var cuboid in cuboids) totalCount += cuboid.ToShaderData().Length;

        if (totalCount == 0)
        {
            material.SetInt(QuadCountPropertyID, 0);
            return;
        }

        var quadData = new ShaderStructs.Quad[totalCount];
        int index = 0;

        // Add all quads
        foreach (var quad in quads)
        {
            var data = quad.ToShaderData();
            for (int i = 0; i < data.Length; i++)
                quadData[index++] = data[i];
        }

        foreach (var plane in planes)
        {
            var data = plane.ToShaderData();
            for (int i = 0; i < data.Length; i++)
                quadData[index++] = data[i];
        }

        foreach (var cuboid in cuboids)
        {
            var data = cuboid.ToShaderData();
            for (int i = 0; i < data.Length; i++)
                quadData[index++] = data[i];
        }

        _quadBuffer = new ComputeBuffer(quadData.Length, System.Runtime.InteropServices.Marshal.SizeOf<ShaderStructs.Quad>());
        _quadBuffer.SetData(quadData);

        material.SetBuffer(QuadBufferPropertyID, _quadBuffer);
        material.SetInt(QuadCountPropertyID, quadData.Length);
    }

    void ReleaseBuffers()
    {
        _sphereBuffer?.Release();
        _quadBuffer?.Release();
    }

    void ReleaseTextures()
    {
        if (_accumulatedTexture != null)
        {
            _accumulatedTexture.Release();
            _accumulatedTexture = null;
        }
    }
}