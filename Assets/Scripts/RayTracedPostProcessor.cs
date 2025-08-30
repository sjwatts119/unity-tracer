// Assets/Scripts/RayTracedPostProcessor.cs
using System;
using Core;
using UnityEngine;
using Shaders;

[ExecuteAlways, ImageEffectAllowedInSceneView]
public class RayTracedPostProcessor : MonoBehaviour
{
    [Header("Info")] public int accumulatedFrames;

    [Header("Anti-Aliasing")]
    [Range(1, 500)]
    public int samplesPerPixel = 1;
    
    [Header("Ray Tracing")]
    [Range(1, 250)]
    public int rayMaxDepth = 10;
    
    [Header("Rendering")]
    public bool accumulateFrames = true;
    
    [Header("Camera")]
    [Range(0.1f, 50f)]
    public float cameraFocalDistance = 1.0f;
    
    [Range(0.0f, 100f)]
    public float cameraDefocusAngle = 0f;

    [Header("Accumulation")]
    [SerializeField] private Shader accumulationShader;

    [Header("Export Render")] 
    [SerializeField] private int exportAtFrame;
    
    // Core components
    private SceneDataProvider _sceneProvider;
    private BufferManager _bufferManager;
    private FrameAccumulator _accumulator;
    private Material _rayMaterial;
    private Camera _camera;
    
    
    // Cached IDs   
    private readonly int _samplesPerPixelID = Shader.PropertyToID("SamplesPerPixel");
    private readonly int _rayMaxDepthID = Shader.PropertyToID("RayMaxDepth");
    private readonly int _cameraFocalDistanceID = Shader.PropertyToID("CameraFocalDistance");
    private readonly int _cameraDefocusAngleID = Shader.PropertyToID("CameraDefocusAngle");
    private readonly int _cameraPlaneWidthID = Shader.PropertyToID("CameraPlaneWidth");
    private readonly int _cameraPlaneHeightID = Shader.PropertyToID("CameraPlaneHeight");
    private readonly int _cameraLocalToWorldID = Shader.PropertyToID("CameraLocalToWorld");
    
    private void Start()
    {
        Initialise();
    }
    
    private void OnEnable()
    {
        Initialise();
    }
    
    private void Initialise()
    {
        _sceneProvider = new SceneDataProvider();
        _bufferManager = new BufferManager();
        _accumulator = new FrameAccumulator(accumulationShader);
        _rayMaterial = RayShader.Material;
        _camera = GetComponent<Camera>();
    }
    
    private void OnValidate()
    {
        if (!Application.isPlaying) return;
        
        _accumulator?.Reset();
    }
    
    private void OnDisable()
    {
        Cleanup();
    }
    
    private void OnDestroy()
    {
        Cleanup();
    }
    
    private void Cleanup()
    {
        _bufferManager?.Dispose();
        _accumulator?.Dispose();
        _bufferManager = null;
        _accumulator = null;
    }
    
    public void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        var sceneData = _sceneProvider.GatherSceneData();
        
        _bufferManager.CreateBuffersFromSceneData(sceneData, _rayMaterial);
        
        // Configure shader parameters
        ConfigureCameraParameters(_rayMaterial);
        ConfigureRayTracingParameters(_rayMaterial);

        if (!accumulateFrames || !Application.isPlaying)
        {
            _accumulator.RenderSingleFrame(source, destination, _rayMaterial);
            accumulatedFrames = _accumulator.FrameCount;
            return;
        }
        
        _accumulator.AccumulateFrame(source, destination, _rayMaterial);
        accumulatedFrames = _accumulator.FrameCount;
        
        if (exportAtFrame > 0 && accumulatedFrames == exportAtFrame)
        {
            ExportRender();
        }
    }
    
    private void ConfigureCameraParameters(Material material)
    {
        // Calculate the dimensions of the camera's near plane
        var planeHeight = 2.0f * cameraFocalDistance * Mathf.Tan(_camera.fieldOfView * 0.5f * Mathf.Deg2Rad);
        var planeWidth = planeHeight * _camera.aspect;
        
        material.SetFloat(_cameraFocalDistanceID, cameraFocalDistance);
        material.SetFloat(_cameraDefocusAngleID, cameraDefocusAngle);
        material.SetFloat(_cameraPlaneWidthID, planeWidth);
        material.SetFloat(_cameraPlaneHeightID, planeHeight);
        material.SetMatrix(_cameraLocalToWorldID, _camera.transform.localToWorldMatrix);
    }
    
    private void ConfigureRayTracingParameters(Material material)
    {
        material.SetInt(_samplesPerPixelID, samplesPerPixel);
        material.SetInt(_rayMaxDepthID, rayMaxDepth);
    }
    
    private static void ExportRender()
    {
        const string directoryName = "Renders";
        var fileName = "Render_" + DateTime.Now.ToString("yyyy-MM-dd_HH-mm-ss") + ".png";
        var fullPath = Application.dataPath + "/" + directoryName + "/" + fileName;
        
        ScreenCapture.CaptureScreenshot(fullPath);
        Debug.Log($"Exported render to: {fullPath}");
    }
}