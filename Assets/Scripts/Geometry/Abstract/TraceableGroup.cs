using Geometry.Interfaces;
using Geometry.Structs;

namespace Geometry.Abstract
{
    public abstract class TraceableGroup : Traceable
    {
        public abstract IPrimitive[] ToPrimitives();

        public abstract PrimitiveGroup ToGroup();
    }
}