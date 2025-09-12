using System;
using System.Collections.Generic;
using System.Linq;
using Geometry.Abstract;
using Geometry.Structs;
using UnityEngine;
using Object = UnityEngine.Object;

namespace Core
{
public class SceneDataProvider
{
    public SceneData GatherSceneData()
    {
        var data = new SceneData();
        
        // Get all individual primitives
        var spheres = Object.FindObjectsByType<Geometry.Sphere>(FindObjectsSortMode.None);
        var quads = Object.FindObjectsByType<Geometry.Quad>(FindObjectsSortMode.None);
        var planes = Object.FindObjectsByType<Geometry.Plane>(FindObjectsSortMode.None);
        var cuboids = Object.FindObjectsByType<Geometry.Cuboid>(FindObjectsSortMode.None);
        var meshes = Object.FindObjectsByType<Geometry.Mesh>(FindObjectsSortMode.None);
        
        // Convert to struct arrays
        var spherePrimitives = spheres.Select(sphere => (Sphere)sphere.ToPrimitive()).ToArray();
        var quadPrimitives = quads.Select(quad => (Quad)quad.ToPrimitive())
                                 .Concat(planes.Select(plane => (Quad)plane.ToPrimitive()))
                                 .ToArray();
        var cuboidPrimitives = cuboids.Select(cuboid => (Cuboid)cuboid.ToPrimitive()).ToArray();
        var trianglePrimitives = meshes.SelectMany(mesh => mesh.GetPrimitives()).ToArray();
        
        data.spheres = spherePrimitives;
        data.quads = quadPrimitives;
        data.cuboids = cuboidPrimitives;
        data.triangles = trianglePrimitives;
        
        return data;
    }
}
}