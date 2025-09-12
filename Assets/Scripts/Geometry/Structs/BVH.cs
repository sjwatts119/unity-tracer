using Unity.Mathematics;

namespace Geometry.Structs
{
    public class Bvh
    {
        public struct BvhNode
        {
            public float3 aabbMin;
            public float3 aabbMax;
            public int leftFirst; // If leaf node, index of first primitive. If internal node, index of first child node.
            public int primitiveCount; // If leaf node, number of primitives. If internal node, 0.
            bool isLeaf => primitiveCount > 0;
        }

        public struct BvhBin
        {
            public AABB aabb;
            public int count;
            
            public static BvhBin CreateEmpty()
            {
                return new BvhBin
                {
                    aabb = AABB.CreateEmpty(),
                    count = 0
                };
            }
        }
    }
}