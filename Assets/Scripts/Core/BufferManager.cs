using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
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
        private static readonly int GroupedSphereCountID = Shader.PropertyToID("GroupedSphereCount");
        private static readonly int GroupedQuadCountID = Shader.PropertyToID("GroupedQuadCount");
        private static readonly int GroupedCuboidCountID = Shader.PropertyToID("GroupedCuboidCount");
        private static readonly int GroupedTriangleCountID = Shader.PropertyToID("GroupedTriangleCount");
        private static readonly int PrimitiveGroupBufferID = Shader.PropertyToID("PrimitiveGroupBuffer");
        private static readonly int PrimitiveGroupCountID = Shader.PropertyToID("PrimitiveGroupCount");
        
        public void CreateBuffersFromSceneData(SceneData data, Material material)
        {
            // Clean up old buffers
            Dispose();
            
            // Create buffers for individual primitives
            CreateBuffer("Spheres", data.individualSpheres.Concat(data.groupedSpheres).ToArray(), 
                material, SphereBufferID, SphereCountID);
            CreateBuffer("Quads", data.individualQuads.Concat(data.groupedQuads).ToArray(), 
                material, QuadBufferID, QuadCountID);
            CreateBuffer("Cuboids", data.individualCuboids.Concat(data.groupedCuboids).ToArray(), 
                material, CuboidBufferID, CuboidCountID);
            CreateBuffer("Triangles", data.individualTriangles.Concat(data.groupedTriangles).ToArray(), 
                material, TriangleBufferID, TriangleCountID);
            
            // Set grouped primitive counts
            material.SetInt(GroupedSphereCountID, data.groupedSpheres.Length);
            material.SetInt(GroupedQuadCountID, data.groupedQuads.Length);
            material.SetInt(GroupedCuboidCountID, data.groupedCuboids.Length);
            material.SetInt(GroupedTriangleCountID, data.groupedTriangles.Length);
            
            // Create primitive group buffer
            CreateBuffer("PrimitiveGroups", data.primitiveGroups, 
                material, PrimitiveGroupBufferID, PrimitiveGroupCountID);
            
            // TODO add BVH node buffers
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