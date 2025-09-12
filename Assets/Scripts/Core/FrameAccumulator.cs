using System;
using UnityEngine;

namespace Core
{
    public class FrameAccumulator : IDisposable
    {
        private RenderTexture _accumulatedTexture;
        private Material _accumulationMaterial;
        private int _frameCount = 0;
        
        private readonly int _frameNumberID = Shader.PropertyToID("FrameNumber");
        private readonly int _previousTextureID = Shader.PropertyToID("PreviousTexture");

        public int FrameCount => _frameCount;

        public FrameAccumulator(Shader accumulationShader)
        {
            if (accumulationShader == null) return;
            
            _accumulationMaterial = new Material(accumulationShader);
        }

        public void Reset()
        {
            _frameCount = 0;
            
            if (_accumulatedTexture == null) return;
            
            RenderTexture.active = _accumulatedTexture;
            GL.Clear(true, true, Color.black);
            RenderTexture.active = null;
        }

        public void AccumulateFrame(
            RenderTexture source,
            RenderTexture destination,
            Material rayMaterial
        ) {
            EnsureTextureExists(source.width, source.height);

            // Copy previous accumulated frame
            var prevFrame = RenderTexture.GetTemporary(source.width, source.height, 0, RenderTextureFormat.ARGBFloat);
            Graphics.Blit(_accumulatedTexture, prevFrame);

            // Render current frame with ray tracing
            rayMaterial.SetInt(_frameNumberID, _frameCount);
            var currentFrame = RenderTexture.GetTemporary(source.width, source.height, 0, RenderTextureFormat.ARGBFloat);
            Graphics.Blit(source, currentFrame, rayMaterial);

            // Blend frames using accumulation shader
            if (_accumulationMaterial != null)
            {
                _accumulationMaterial.SetInt(_frameNumberID, _frameCount);
                _accumulationMaterial.SetTexture(_previousTextureID, prevFrame);
                Graphics.Blit(currentFrame, _accumulatedTexture, _accumulationMaterial);
            }
            else
            {
                // If no accumulation shader, blit the current frame directly
                Graphics.Blit(currentFrame, _accumulatedTexture);
            }

            // Output accumulated result to screen
            Graphics.Blit(_accumulatedTexture, destination);

            // Cleanup temporary textures
            RenderTexture.ReleaseTemporary(currentFrame);
            RenderTexture.ReleaseTemporary(prevFrame);

            _frameCount++;
        }

        public void RenderSingleFrame(RenderTexture source, RenderTexture destination, Material rayMaterial)
        {
            Reset();
            rayMaterial.SetInt(_frameNumberID, 0);
            Graphics.Blit(source, destination, rayMaterial);
        }

        private void EnsureTextureExists(int width, int height)
        {
            // If the texture already exists and matches the required size, do nothing
            if (_accumulatedTexture != null && _accumulatedTexture.width == width && _accumulatedTexture.height == height) return;

            _accumulatedTexture?.Release();
            _accumulatedTexture = new RenderTexture(width, height, 0, RenderTextureFormat.ARGBFloat)
            {
                enableRandomWrite = true,
            };
            _accumulatedTexture.Create();
            Reset();
        }

        public void Dispose()
        {
            _accumulatedTexture?.Release();
            _accumulatedTexture = null;

            if (_accumulationMaterial == null) return;
            
            UnityEngine.Object.DestroyImmediate(_accumulationMaterial);
            _accumulationMaterial = null;
        }
    }
}