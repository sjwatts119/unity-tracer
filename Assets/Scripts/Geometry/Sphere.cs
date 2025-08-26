using Geometry.Abstract;
using Geometry.Interfaces;
using UnityEngine;
using Vector3 = UnityEngine.Vector3;

namespace Geometry
{
    public class Sphere : TraceablePrimitive
    {
        private Vector3 LocalScale => transform.localScale;
        
        private Vector3 LocalPosition => transform.localPosition;
        
        private float ScaledRadius => Mathf.Max(LocalScale.x, LocalScale.y, LocalScale.z) * 0.5f;
        
        public override IPrimitive ToPrimitive()
        {
            return new Geometry.Structs.Sphere
            {
                centre = LocalPosition,
                radius = ScaledRadius,
                material = GetMaterial()
            };
        }
    }
}