using Geometry.Abstract;
using Geometry.Interfaces;
using UnityEngine;
using Vector3 = UnityEngine.Vector3;

namespace Geometry
{
    public class Sphere : TraceablePrimitive
    {
        private Vector3 WorldScale => transform.lossyScale;
    
        private Vector3 WorldPosition => transform.position;
    
        private float ScaledRadius => Mathf.Max(WorldScale.x, WorldScale.y, WorldScale.z) * 0.5f;
    
        public override IPrimitive ToPrimitive()
        {
            return new Geometry.Structs.Sphere
            {
                centre = WorldPosition,
                radius = ScaledRadius,
                material = GetMaterial()
            };
        }
    }
}