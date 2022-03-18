import numpy as np
from foronoi import Voronoi, Polygon, Visualizer, VoronoiObserver
import pandas as pd

obs_data = pd.read_excel('C:\JunXiang\傑明工程\Fem\excel\\Obs_conc_v1.xlsx')
def voronoi(obs_data):
    """
    Ref : https://github.com/Yatoom/foronoi
    """
    area = {}
    points = []
    for i in range(len(obs_data)):
        points.append([obs_data['Y'][i], obs_data['Z'][i]])

    points = np.array(points)

    # Boundary
    polygon = Polygon([
    (0, 0),
    (50, 0),
    (0, 20),
    (50, 20)
    ])

    v = Voronoi(polygon)

    v.attach_observer(VoronoiObserver())
    v.create_diagram(points=points)
    edges = v.edges
    vertices = v.vertices
    arcs = v.arcs

    Visualizer(v, canvas_offset=1).plot_sites(show_labels=False, color='r').plot_edges(show_labels=False).show()

    for i, point in enumerate(v.sites):
        print(f"{point.xy} \t {point.area()}")
        area[obs_data['Node'][i]] = round(point.area(), 1)

    return area
