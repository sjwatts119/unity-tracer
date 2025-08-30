using Geometry.Interfaces;
using Geometry.Structs;

namespace Geometry.Abstract
{
    public abstract class TraceableGroup : Traceable
    {
        public abstract AABB ToAABB();
    
        // Return all primitives in this group, organized by type
        public abstract PrimitiveCollection GetPrimitives();
    }

}