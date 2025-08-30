namespace Geometry.Structs
{
    [System.Serializable]
    public struct PrimitiveCollection
    {
        public Sphere[] spheres;
        public Quad[] quads;
        public Cuboid[] cuboids;
        public Triangle[] triangles;
    }
}