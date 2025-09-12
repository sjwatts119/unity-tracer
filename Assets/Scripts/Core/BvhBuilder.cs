using System;
using System.Collections.Generic;
using System.Diagnostics;
using Geometry.Structs;
using Unity.Mathematics;
using Debug = UnityEngine.Debug;

namespace Core
{
    public class BvhBuilder
    {
        private const int Intervals = 16; // Number of bins to use for SAH evaluation
        
        private readonly IReadOnlyList<Triangle> _triangles;
        private int TriangleCount => _triangles.Count;

        // Pad the triangle indices array to avoid modifying the input triangles list
        private readonly int[] _triangleIndices;

        // We can only have a maximum of 2N-1 nodes in a BVH for N triangles
        private readonly Bvh.BvhNode[] _nodes;
        private const int RootNodeIndex = 0;
        private int _totalNodes = 1;
        
        public BvhBuilder(IReadOnlyList<Triangle> triangles)
        {
            _triangles = triangles;
            
            _triangleIndices = new int[TriangleCount];
            for (var i = 0; i < TriangleCount; i++) _triangleIndices[i] = i;
            
            _nodes = new Bvh.BvhNode[math.max(TriangleCount * 2 - 1, 0)];
        }

        public (Bvh.BvhNode[] nodes, Triangle[] orderedTriangles) Build()
        {
            if (TriangleCount == 0) return (_nodes, Array.Empty<Triangle>());

            Stopwatch stopwatch = Stopwatch.StartNew();
    
            var root = _nodes[RootNodeIndex];
            root.leftFirst = 0;
            root.primitiveCount = TriangleCount;
            _nodes[RootNodeIndex] = root;

            UpdateNodeBounds(RootNodeIndex);
            Subdivide(RootNodeIndex);
    
            stopwatch.Stop();
            Debug.Log($"Built BVH in {stopwatch.ElapsedMilliseconds} ms");

            return (_nodes, GetOrderedTriangles());
        }
        
        private Triangle[] GetOrderedTriangles()
        {
            var ordered = new Triangle[TriangleCount];
            for (var i = 0; i < TriangleCount; i++) ordered[i] = _triangles[_triangleIndices[i]];
            return ordered;
        }

        private void UpdateNodeBounds(int nodeIndex)
        {
            Bvh.BvhNode node = _nodes[nodeIndex];

            // Set the first AABB bounds to max and min values to start
            node.aabbMin = new float3(float.MaxValue);
            node.aabbMax = new float3(float.MinValue);

            // Iterate over all triangles in this node and expand the AABB to include them
            for (int first = node.leftFirst, i = 0; i < node.primitiveCount; i++)
            {
                var leafTriangle = _triangles[_triangleIndices[first + i]];
                
                // Expand the AABB to include the triangle's vertices
                // We should refactor this into a utility function on an AABB struct really but this is fine for now
                node.aabbMin = math.min(node.aabbMin, leafTriangle.v0);
                node.aabbMin = math.min(node.aabbMin, leafTriangle.v1);
                node.aabbMin = math.min(node.aabbMin, leafTriangle.v2);
                node.aabbMax = math.max(node.aabbMax, leafTriangle.v0);
                node.aabbMax = math.max(node.aabbMax, leafTriangle.v1);
                node.aabbMax = math.max(node.aabbMax, leafTriangle.v2);
            }
            
            // Write back the updated node
            _nodes[nodeIndex] = node;
        }

        private (float bestCost, int axis, float splitPosition) FindBestSplitPlane(Bvh.BvhNode node)
        {
            var bestCost = float.MaxValue;
            var bestAxis = 0;
            var bestSplitPosition = 0f;
            
            // Iterate over each axis
            for(var x = 0; x < 3; x++)
            {
                float min = float.MaxValue;
                float max = float.MinValue;
                
                // Find the min and max centroid along this axis
                for (var y = 0; y < node.primitiveCount; y++)
                {
                    var triangle = _triangles[_triangleIndices[node.leftFirst + y]];
                    min = math.min(min, triangle.centroid[x]);
                    max = math.max(max, triangle.centroid[x]);
                }

                if (math.abs(max - min) < 1e-5) continue; // No point splitting along this axis if the extent is zero
                
                Bvh.BvhBin[] bins = new Bvh.BvhBin[Intervals];
                for (int b = 0; b < Intervals; b++) bins[b] = Bvh.BvhBin.CreateEmpty();
                
                var scale = Intervals / (max - min);
                
                // Build up the bins with triangle counts and AABBs
                for (var i = 0; i < node.primitiveCount; i++)
                {
                    var triangle = _triangles[_triangleIndices[node.leftFirst + i]];
                    var binIndex = math.min(Intervals - 1, (int)((triangle.centroid[x] - min) * scale));
                    bins[binIndex].count++;
                    bins[binIndex].aabb.ExpandToInclude(triangle.v0);
                    bins[binIndex].aabb.ExpandToInclude(triangle.v1);
                    bins[binIndex].aabb.ExpandToInclude(triangle.v2);
                }
                
                float[] leftArea = new float[Intervals - 1];
                float[] rightArea = new float[Intervals - 1];
                int[] leftCount = new int[Intervals - 1];
                int[] rightCount = new int[Intervals - 1];
                AABB leftAABB = AABB.CreateEmpty();
                AABB rightAABB = AABB.CreateEmpty();
                var leftSum = 0;
                var rightSum = 0;

                for (var i = 0; i < Intervals - 1; i++)
                {
                    leftSum += bins[i].count;
                    leftCount[i] = leftSum;
                    leftAABB.ExpandToInclude(bins[i].aabb);
                    leftArea[i] = leftAABB.SurfaceArea();
                    
                    rightSum += bins[Intervals - 1 - i].count;
                    rightCount[Intervals - 2 - i] = rightSum;
                    rightAABB.ExpandToInclude(bins[Intervals - 1 - i].aabb);
                    rightArea[Intervals - 2 - i] = rightAABB.SurfaceArea();
                }

                scale = (max - min) / Intervals;
                
                // Evaluate the cost for each split position
                for (var i = 0; i < Intervals - 1; i++)
                {
                    var cost = leftCount[i] * leftArea[i] + rightCount[i] * rightArea[i];
                    
                    if (cost >= bestCost) continue;
                    
                    bestAxis = x;
                    bestSplitPosition = min + scale * (i + 1);
                    bestCost = cost;
                }
            }

            return (bestCost, bestAxis, bestSplitPosition);
        }

        private void Subdivide(int nodeIndex)
        {
            // Get the node to subdivide
            var node = _nodes[nodeIndex];

            // Determine the axis and position to split along
            var (bestCost, axis, splitPos) = FindBestSplitPlane(node);
            
            var parentCost = EvaluateLeafCost(node);
            
            if (bestCost >= parentCost) return; // Best split found has a higher cost than not splitting, so don't split

            // Partition the triangles into two groups based on the split position
            var i = node.leftFirst;
            var j = i + node.primitiveCount - 1;
            while (i <= j)
            {
                var triangle = _triangles[_triangleIndices[i]];
                if (triangle.centroid[axis] < splitPos)
                {
                    i++;
                }
                else
                {
                    // Swap the triangles to the other side
                    (_triangleIndices[i], _triangleIndices[j]) = (_triangleIndices[j--], _triangleIndices[i]);
                }
            }

            var leftCount = i - node.leftFirst;
            
            // If one of the sides is empty, we can't split this node, just return this as a leaf
            if (leftCount == 0 || leftCount == node.primitiveCount) return;
            
            // Create our new child nodes
            var leftChildIndex = _totalNodes++;
            var rightChildIndex = _totalNodes++;
            
            // Set up the left child node
            _nodes[leftChildIndex] = new Bvh.BvhNode
            {
                leftFirst = node.leftFirst,
                primitiveCount = leftCount
            };
            
            // Set up the right child node
            _nodes[rightChildIndex] = new Bvh.BvhNode
            {
                leftFirst = i,
                primitiveCount = node.primitiveCount - leftCount
            };
            
            // Update the parent node to point to the children and to be a non-leaf node
            node.leftFirst = leftChildIndex;
            node.primitiveCount = 0;
            
            // Write back the updated parent node
            _nodes[nodeIndex] = node;
            
            // Recursively update AABB bounds for the child nodes
            UpdateNodeBounds(leftChildIndex);
            UpdateNodeBounds(rightChildIndex);
            
            // Recursively subdivide the child nodes
            Subdivide(leftChildIndex);
            Subdivide(rightChildIndex);
        }
        
        /**
         * Evaluate the cost of not splitting a node (i.e. making it a leaf).
         */
        float EvaluateLeafCost(Bvh.BvhNode node)
        {
            var extents = node.aabbMax - node.aabbMin;
            var surfaceArea = extents.x * extents.y + extents.y * extents.z + extents.z * extents.x;
            return node.primitiveCount * surfaceArea;
        }
    }
}