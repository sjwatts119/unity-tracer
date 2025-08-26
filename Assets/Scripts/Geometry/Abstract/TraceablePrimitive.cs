using Geometry.Interfaces;

namespace Geometry.Abstract
{
    public abstract class TraceablePrimitive : Traceable
    {
        public abstract IPrimitive ToPrimitive();
    }
}