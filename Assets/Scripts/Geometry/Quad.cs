using Geometry.Abstract;
using Geometry.Interfaces;
using Vector3 = UnityEngine.Vector3;

namespace Geometry
{
    public class Quad : TraceablePrimitive
    {
        private Vector3 ScaledRight => transform.right * transform.localScale.x;
        
        private Vector3 ScaledUp => transform.up * transform.localScale.y;
        
        private Vector3 Corner => transform.position - (ScaledRight * 0.5f) - (ScaledUp * 0.5f);
        
        public override IPrimitive ToPrimitive()
        {
            return new Geometry.Structs.Quad
            {
                q = Corner,
                u = ScaledRight,
                v = ScaledUp,
                material = GetMaterial()
            };
        }
    }
}
