# pass facets from the normal tbem structures to this module
# receiving facets, process them into a mesh
# find extents of box from the facets
import meshpy.triangle as triangle
import numpy as np

def determine_extents(facets):
    min_corner = np.min(np.min(facets, axis = 0), axis = 0)
    max_corner = np.max(np.max(facets, axis = 0), axis = 0)
    return min_corner, max_corner

def expand_extents(min_corner, max_corner, factor = 2.0):
    center = (min_corner + max_corner) / 2.0
    dist = (max_corner - min_corner) / 2.0
    return center - dist * factor, center + dist * factor

def extents_to_box_2d(min_corner, max_corner):
    pts = [
        [min_corner[0], min_corner[1]],
        [max_corner[0], min_corner[1]],
        [max_corner[0], max_corner[1]],
        [min_corner[0], max_corner[1]]
    ]
    return np.array([
        [pts[0], pts[1]],
        [pts[1], pts[2]],
        [pts[2], pts[3]],
        [pts[3], pts[0]]
    ])

def add_extent_surface(input_facets):
    extents_box = extents_to_box_2d(
        *expand_extents(*determine_extents(input_facets))
    )
    return np.concatenate((input_facets, extents_box))

def create_mesh(facets):
    meshpy_vs = facets.reshape((facets.shape[0] * facets.shape[1], facets.shape[2]))
    meshpy_facets = np.arange(meshpy_vs.shape[0]).reshape(
        (facets.shape[0], facets.shape[1])
    )
    meshpy_markers = np.arange(meshpy_facets.shape[0])
    info = triangle.MeshInfo()
    info.set_points(meshpy_vs)
    info.set_facets(meshpy_facets, facet_markers = meshpy_markers)
    mesh = triangle.build(info, generate_faces = True)
    pts = np.array(mesh.points)
    tris = np.array(mesh.elements)
    import ipdb; ipdb.set_trace()
    import matplotlib.pyplot as plt
    plt.triplot(pts[:, 0], pts[:, 1], tris)
    plt.show()

