using System.IO;
using UnityEngine;

namespace Core
{
    public class RayShader
    {
        public static readonly string Name = "RayTracer/RayShader";

        public static Material GetMaterial()
        {
            return new Material(Shader.Find(Name));
        }
    }
}