using Geometry.Abstract;
using Geometry.Interfaces;
using UnityEngine;

namespace Geometry
{
    public class Cuboid : TraceablePrimitive
    {
        // Get the rotation matrix from the transform's rotation
        private Matrix4x4 RotationMatrix => Matrix4x4.Rotate(transform.rotation);
        
        // Create worldToLocal and localToWorld matrices from the rotation matrix
        private Matrix4x4 WorldToLocal => RotationMatrix.transpose;
        
        private Matrix4x4 LocalToWorld => RotationMatrix;
        
        public override IPrimitive ToPrimitive()
        {
            return new Geometry.Structs.Cuboid
            {
                centre = transform.position,
                size = transform.localScale,
                worldToLocal = WorldToLocal,
                localToWorld = LocalToWorld,
                material = GetMaterial()
            };
        }
    }
}