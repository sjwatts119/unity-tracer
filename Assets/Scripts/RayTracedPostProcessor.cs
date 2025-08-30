using System;
using System.Linq;
using System.Runtime.InteropServices;
using UnityEngine;
using Geometry.Abstract;
using Geometry.Interfaces;
using Geometry.Structs;
using Shaders;
using Unity.VisualScripting;
using Cuboid = Geometry.Cuboid;
using Plane = Geometry.Plane;
using Quad = Geometry.Quad;
using Sphere = Geometry.Sphere;

[ExecuteAlways, ImageEffectAllowedInSceneView]
public class RayTracedPostProcessor : MonoBehaviour
{
    [Header("Info")] 
    public int renderedFrames;
    
    [Header("Anti-Aliasing")]
    [Range(1, 500)]
    public int samplesPerPixel;
    
    [Header("Ray Tracing")]
    [Range(1, 250)]
    public int rayMaxDepth;
    
    [Header("Rendering")]
    public bool accumulateFrames = true;
    
    [Header("Camera")]
    [Range(0.1f, 50f)]
    public float cameraFocalDistance = 1.0f;
    
    [Range(0.0f, 100f)]
    public float cameraDefocusAngle = 10f;

    [Header("Accumulation")]
    [SerializeField] private Shader accumulationShader;

    [Header("Export Render")] 
    [SerializeField] private int exportAtFrame;
    
    // Geometric data buffers
    private ComputeBuffer _sphereBuffer;
    private ComputeBuffer _quadBuffer;
    private ComputeBuffer _cuboidBuffer;
    private ComputeBuffer _triangleBuffer;
    private ComputeBuffer _groupBuffer;
    
    // Frame accumulation data
    private RenderTexture _accumulatedTexture;
    private Material _accumulationMaterial;
    
    // Primitive Buffers
    private static readonly int SphereBufferPropertyID = Shader.PropertyToID("SphereBuffer");
    private static readonly int QuadBufferPropertyID = Shader.PropertyToID("QuadBuffer");
    private static readonly int CuboidBufferPropertyID = Shader.PropertyToID("CuboidBuffer");
    private static readonly int TriangleBufferPropertyID = Shader.PropertyToID("TriangleBuffer");
    
    // Primitive Counts
    private static readonly int SphereCountPropertyID = Shader.PropertyToID("SphereCount");
    private static readonly int QuadCountPropertyID = Shader.PropertyToID("QuadCount");
    private static readonly int CuboidCountPropertyID = Shader.PropertyToID("CuboidCount");
    private static readonly int TriangleCountPropertyID = Shader.PropertyToID("TriangleCount");
    
    // Grouped Primitive Counts
    private static readonly int GroupedSphereCountPropertyID = Shader.PropertyToID("GroupedSphereCount");
    private static readonly int GroupedQuadCountPropertyID = Shader.PropertyToID("GroupedQuadCount");
    private static readonly int GroupedCuboidCountPropertyID = Shader.PropertyToID("GroupedCuboidCount");
    private static readonly int GroupedTriangleCountPropertyID = Shader.PropertyToID("GroupedTriangleCount");

    // Primitive Groups
    private static readonly int PrimitiveGroupBufferPropertyID = Shader.PropertyToID("PrimitiveGroupBuffer");
    private static readonly int PrimitiveGroupCountPropertyID = Shader.PropertyToID("PrimitiveGroupCount");
    
    // Camera
    private static readonly int CameraFocalDistancePropertyID = Shader.PropertyToID("CameraFocalDistance");
    private static readonly int CameraPlaneWidthPropertyID = Shader.PropertyToID("CameraPlaneWidth");
    private static readonly int CameraPlaneHeightPropertyID = Shader.PropertyToID("CameraPlaneHeight");
    private static readonly int CameraDefocusAnglePropertyID = Shader.PropertyToID("CameraDefocusAngle");
    private static readonly int CameraLocalToWorldPropertyID = Shader.PropertyToID("CameraLocalToWorld");
    
    // Anti-aliasing
    private static readonly int SamplesPerPixelPropertyID = Shader.PropertyToID("SamplesPerPixel");
    private static readonly int RayMaxDepthPropertyID = Shader.PropertyToID("RayMaxDepth");
    
    // Frame accumulation
    private static readonly int FrameNumberPropertyID = Shader.PropertyToID("FrameNumber");
    private static readonly int PreviousTexturePropertyID = Shader.PropertyToID("PreviousTexture");

    private void Start()
    {
        renderedFrames = 0;
        InitAccumulation();
    }
    
    void OnDisable()
    {
        ReleaseBuffers();
        ReleaseAccumulatedTextures();
    }

    void OnDestroy()
    {
        ReleaseBuffers();
        ReleaseAccumulatedTextures();
    }
    
    // Reset accumulation when parameters change
    void OnValidate()
    {
        if (!Application.isPlaying) return;
        
        renderedFrames = 0;
        ReleaseAccumulatedTextures();
    }

    // Called after all rendering is complete to render image
    public void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        var rayMaterial = RayShader.Material;
        
        PopulateAntialiasingData(rayMaterial);
        PopulateRaytracingData(rayMaterial);
        PopulateCameraData(rayMaterial);
        PopulateRenderingData(rayMaterial);
        
        PopulateBufferData<Sphere, Geometry.Structs.Sphere>(ref _sphereBuffer, rayMaterial, SphereBufferPropertyID, SphereCountPropertyID);
        PopulateBufferData<Cuboid, Geometry.Structs.Cuboid>(ref _cuboidBuffer, rayMaterial, CuboidBufferPropertyID, CuboidCountPropertyID);
        
        PopulateQuadData(rayMaterial);
        
        // We should populate groups last as it relies on indexes from other buffers
        PopulateMeshData(rayMaterial);
        
        // If we are not accumulating frames, just render the current frame
        if (!accumulateFrames)
        {
            renderedFrames = 0;
            DrawRawFrame(source, destination, rayMaterial);
            return;
        }

        InitFrame(source);
        
        DrawAccumulatedFrame(source, destination, rayMaterial);

        if (exportAtFrame > 0 && renderedFrames == exportAtFrame)
        {
            ExportRender();
        }

        renderedFrames += Application.isPlaying ? 1 : 0;
    }

    void InitFrame(RenderTexture source)
    {
        InitAccumulation();
        CreateAccumulationTexture(source.width, source.height);
    }

    private void InitAccumulation()
    {
        // Create accumulation material if needed
        if (accumulationShader == null || _accumulationMaterial != null) return;
    
        _accumulationMaterial = new Material(accumulationShader);
    }

    private void CreateAccumulationTexture(int width, int height)
    {
        // Create or recreate the accumulation texture if size has changed
        if (_accumulatedTexture != null && _accumulatedTexture.width == width && _accumulatedTexture.height == height) return;
    
        if (_accumulatedTexture != null)
        {
            _accumulatedTexture.Release();
        }

        _accumulatedTexture = new RenderTexture(width, height, 0, RenderTextureFormat.ARGBFloat);
        _accumulatedTexture.enableRandomWrite = true;
        _accumulatedTexture.Create();
    }
    
    void DrawRawFrame(RenderTexture source, RenderTexture destination, Material rayMaterial)
    {
        Graphics.Blit(source, destination, rayMaterial);
    }
    
    void DrawAccumulatedFrame(RenderTexture source, RenderTexture destination, Material rayMaterial)
    {
        // Create temporary copy of the previous accumulated frame
        RenderTexture prevFrameCopy = RenderTexture.GetTemporary(source.width, source.height, 0, RenderTextureFormat.ARGBFloat);
        Graphics.Blit(_accumulatedTexture, prevFrameCopy);

        // Render the current frame with ray tracing
        rayMaterial.SetInt(FrameNumberPropertyID, renderedFrames);
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

        // Draw the accumulated result to screen
        Graphics.Blit(_accumulatedTexture, destination);

        // Clean up temporary textures
        RenderTexture.ReleaseTemporary(currentFrame);
        RenderTexture.ReleaseTemporary(prevFrameCopy);
    }
    
    /*
     * Exporting renders
     */

    private void ExportRender()
    {
        string directoryName = "Renders";
        string fileName = "Render_" + DateTime.Now.ToString("yyyy-MM-dd_HH-mm-ss") + ".png";
        string fullPath = Application.dataPath + "/" + directoryName + "/" + fileName;
        
        ScreenCapture.CaptureScreenshot(fullPath);
    }

    /*
     * Populating shader data
     */
    
    private void PopulateCameraData(Material material)
    {
        var cam = GetComponent<Camera>();
        
        // Calculate the dimensions of the camera's near plane
        var planeHeight = 2.0f * cameraFocalDistance * Mathf.Tan(cam.fieldOfView * 0.5f * Mathf.Deg2Rad);
        var planeWidth = planeHeight * cam.aspect;
        
        // Send camera parameters to shader
        material.SetFloat(CameraFocalDistancePropertyID, cameraFocalDistance);
        material.SetFloat(CameraDefocusAnglePropertyID, cameraDefocusAngle);
        material.SetFloat(CameraPlaneWidthPropertyID, planeWidth);
        material.SetFloat(CameraPlaneHeightPropertyID, planeHeight);
        material.SetMatrix(CameraLocalToWorldPropertyID, cam.transform.localToWorldMatrix);
    }
    
    private void PopulateAntialiasingData(Material material)
    {
        material.SetInt(SamplesPerPixelPropertyID, samplesPerPixel);
    }
    
    private void PopulateRaytracingData(Material material)
    {
        material.SetInt(RayMaxDepthPropertyID, rayMaxDepth);
    }
    
    private void PopulateRenderingData(Material material)
    {
        material.SetInt(FrameNumberPropertyID, accumulateFrames ? renderedFrames : 0);
    }

    // Generic method to populate buffer data for different primitive types
    private static void PopulateBufferData<TComponent, TPrimitive>(
        ref ComputeBuffer buffer,
        Material material,
        int bufferPropertyID,
        int countPropertyID)
        where TComponent : TraceablePrimitive
        where TPrimitive : struct, IPrimitive
    {
        var components = FindObjectsByType<TComponent>(FindObjectsSortMode.None);
    
        // Clean up existing buffer
        buffer?.Release();
        buffer = null;
    
        if(components.Length == 0) 
        {
            material.SetInt(countPropertyID, 0);
            return;
        }
    
        var primitiveData = new TPrimitive[components.Length];
        for (var i = 0; i < components.Length; i++)
        {
            primitiveData[i] = (TPrimitive)components[i].ToPrimitive();
        }

        // Create and populate compute buffer
        buffer = new ComputeBuffer(primitiveData.Length, Marshal.SizeOf<TPrimitive>());
        buffer.SetData(primitiveData);
    
        // Bind buffer to shader
        material.SetBuffer(bufferPropertyID, buffer);
        material.SetInt(countPropertyID, primitiveData.Length);
    }

    // Merge quads and planes into a single buffer as they are treated identically inside the shader
    void PopulateQuadData(Material material)
    {
        var quads = FindObjectsByType<Quad>(FindObjectsSortMode.None);
        var planes = FindObjectsByType<Plane>(FindObjectsSortMode.None);

        _quadBuffer?.Release();
        _quadBuffer = null;

        var totalCount = quads.Length + planes.Length;
   
        if (totalCount == 0)
        {
            material.SetInt(QuadCountPropertyID, 0);
            return;
        }

        var quadData = new Geometry.Structs.Quad[totalCount];
        var index = 0;

        // Add all quads
        foreach (var quad in quads)
        {
            quadData[index++] = (Geometry.Structs.Quad)quad.ToPrimitive();
        }

        // Add all planes (converted to quads)
        foreach (var plane in planes)
        {
            quadData[index++] = (Geometry.Structs.Quad)plane.ToPrimitive();
        }

        _quadBuffer = new ComputeBuffer(quadData.Length, Marshal.SizeOf<Geometry.Structs.Quad>());
        _quadBuffer.SetData(quadData);

        material.SetBuffer(QuadBufferPropertyID, _quadBuffer);
        material.SetInt(QuadCountPropertyID, quadData.Length);
    }

    void PopulateMeshData(Material material)
    {
        var meshes = FindObjectsByType<Geometry.Mesh>(FindObjectsSortMode.None);
    
        _triangleBuffer?.Release();
        _groupBuffer?.Release();
        _triangleBuffer = null;
        _groupBuffer = null;
    
        if (meshes.Length == 0)
        {
            material.SetInt(GroupedTriangleCountPropertyID, 0);
            material.SetInt(PrimitiveGroupCountPropertyID, 0);
            return;
        }

        var triangles = new System.Collections.Generic.List<Triangle>();
        var groups = new PrimitiveGroup[meshes.Length];

        foreach (var (mesh, i) in meshes.Select((value, i) => (value, i)))
        {
            // Fix: use current triangle count as start index
            groups[i] = mesh.ToPrimitiveGroup(0, 0, 0, triangles.Count);
        
            triangles.AddRange(mesh.Triangles);
        }
    
        _triangleBuffer = new ComputeBuffer(triangles.Count, Marshal.SizeOf<Triangle>());
        _triangleBuffer.SetData(triangles.ToArray());
        material.SetBuffer(TriangleBufferPropertyID, _triangleBuffer);
        material.SetInt(GroupedTriangleCountPropertyID, triangles.Count);
    
        _groupBuffer = new ComputeBuffer(groups.Length, Marshal.SizeOf<PrimitiveGroup>());
        _groupBuffer.SetData(groups);
        material.SetBuffer(PrimitiveGroupBufferPropertyID, _groupBuffer);
        material.SetInt(PrimitiveGroupCountPropertyID, groups.Length);
    }

    void ReleaseBuffers()
    {
        _sphereBuffer?.Release();
        _quadBuffer?.Release();
        _cuboidBuffer?.Release();
        _triangleBuffer?.Release();
        _groupBuffer?.Release();
    }

    void ReleaseAccumulatedTextures()
    {
        if (_accumulatedTexture == null) return;
        
        _accumulatedTexture.Release();
        _accumulatedTexture = null;
    }
}