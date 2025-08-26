using Geometry.Abstract;
using Geometry.Interfaces;
using Vector3 = UnityEngine.Vector3;

namespace Geometry
{
    public class Plane : TraceablePrimitive
    {
        // Unity's planes are 10x10 versus quads being 1x1, so we just hack in a scale factor lol
        private const float ScaleFactor = 10f;
        
        private Vector3 ScaledRight => transform.right * transform.localScale.x * ScaleFactor;
        
        private Vector3 ScaledForward => transform.forward * transform.localScale.z * ScaleFactor;

        private Vector3 Corner => transform.position - (ScaledRight * 0.5f) - (ScaledForward * 0.5f);

        public override IPrimitive ToPrimitive()
        {
            return new Geometry.Structs.Quad
            {
                q = Corner,
                u = ScaledRight,
                v = ScaledForward,
                material = GetMaterial()
            };
        }
    }
}
