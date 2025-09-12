using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using Geometry.Structs;
using UnityEngine;

namespace Core
{
    public class BufferManager : IDisposable
    {
        private readonly Dictionary<string, ComputeBuffer> _buffers = new();
        
        // Cached shader property IDs
        private static readonly int SphereBufferID = Shader.PropertyToID("SphereBuffer");
        private static readonly int SphereCountID = Shader.PropertyToID("SphereCount");
        private static readonly int QuadBufferID = Shader.PropertyToID("QuadBuffer");
        private static readonly int QuadCountID = Shader.PropertyToID("QuadCount");
        private static readonly int CuboidBufferID = Shader.PropertyToID("CuboidBuffer");
        private static readonly int CuboidCountID = Shader.PropertyToID("CuboidCount");
        private static readonly int TriangleBufferID = Shader.PropertyToID("TriangleBuffer");
        private static readonly int TriangleCountID = Shader.PropertyToID("TriangleCount");
        private static readonly int BvhNodeBufferID = Shader.PropertyToID("BvhNodeBuffer");
        private static readonly int BvhNodeCountID = Shader.PropertyToID("BvhNodeCount");
        
        private (Bvh.BvhNode[] nodes, Triangle[] orderedTriangles) _cachedBvhData;
        
        public void CreateBuffersFromSceneData(SceneData data, Material material)
        {
            // Clean up old buffers
            Dispose();
    
            // If the cache hasn't been built yet, or if the triangles have changed, rebuild the BVH
            if (_cachedBvhData.orderedTriangles == null)
            {
                var bvhBuilder = new BvhBuilder(data.triangles);
                _cachedBvhData = bvhBuilder.Build();
            }
    
            // Create compute buffers for each geometry type
            CreateBuffer("Spheres", data.spheres.ToArray(), 
                material, SphereBufferID, SphereCountID);
            CreateBuffer("Quads", data.quads.ToArray(), 
                material, QuadBufferID, QuadCountID);
            CreateBuffer("Cuboids", data.cuboids.ToArray(), 
                material, CuboidBufferID, CuboidCountID);
    
            // Use the reordered triangles from BVH builder instead of original triangles
            CreateBuffer("Triangles", _cachedBvhData.orderedTriangles, 
                material, TriangleBufferID, TriangleCountID);
            CreateBuffer("BvhNodes", _cachedBvhData.nodes, 
                material, BvhNodeBufferID, BvhNodeCountID);
        }
        
        private void CreateBuffer<T>(string key, T[] data, Material material, 
            int bufferPropertyID, int countPropertyID) where T : struct
        {
            if (data == null || data.Length == 0)
            {
                material.SetInt(countPropertyID, 0);
                return;
            }
            
            var buffer = new ComputeBuffer(data.Length, Marshal.SizeOf<T>());
            buffer.SetData(data);
            _buffers[key] = buffer;
            
            material.SetBuffer(bufferPropertyID, buffer);
            material.SetInt(countPropertyID, data.Length);
        }
        
        public void Dispose()
        {
            foreach (var buffer in _buffers.Values)
            {
                buffer?.Release();
            }
            _buffers.Clear();
        }
    }
}