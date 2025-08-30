using Geometry.Interfaces;
using Geometry.Structs;

namespace Geometry.Abstract
{
    public abstract class TraceableGroup : Traceable
    {
        public abstract AABB ToAABB();
        
        public abstract PrimitiveGroup ToPrimitiveGroup(int currentTris, int currentSpheres, int currentQuads, int currentCuboids);
    }
}