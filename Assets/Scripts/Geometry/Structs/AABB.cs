using Unity.Mathematics;

namespace Geometry.Structs
{
    public struct AABB
    {
        public float3 min;
        public float3 max;

        public AABB(float3 minPos, float3 maxPos)
        {
            min = minPos;
            max = maxPos;
        }

        public static AABB CreateEmpty()
        {
            return new AABB(new float3(float.MaxValue), new float3(float.MinValue));
        }

        public void ExpandToInclude(float3 point)
        {
            min = math.min(min, point);
            max = math.max(max, point);
        }

        public void ExpandToInclude(AABB other)
        {
            min = math.min(min, other.min);
            max = math.max(max, other.max);
        }

        public float SurfaceArea()
        {
            float3 extents = max - min;
            return extents.x * extents.y + extents.y * extents.z + extents.z * extents.x;
        }
    }
}