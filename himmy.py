import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull


#K.P: Parameter initialization 
number_of_points: int = 15
alpha: int = 2
point_set= []
viewpoint_x= float(input("Enter a viewpoint X value:"))
viewpoint_y= float(input("Enter a viewpoint Y value:"))
dimensions: int = 2
reflected_dict = {}
invisible_thickness = 10
# is_3D = (input("Are you using 3 dimensions? (Y/n)"))
# is_3D = is_3D == "Y" or is_3D == "y"
# if is_3D:
#     viewpoint_z = float(input("Enter a viewpoint Z value:"))
#     viewpoint = [viewpoint_x, viewpoint_y, viewpoint_z]
#     dimensions = 3
# else:
#     viewpoint = [viewpoint_x, viewpoint_y]
viewpoint = [viewpoint_x, viewpoint_y]

#K.P: Generate point set
xs = np.linspace(0, 10*np.pi, number_of_points)
for i in xs:
    #TODO: Algo for R3 needed 
    #K.P: This is the circle equation
    # point = [4 + np.cos(i), np.sin(i)]
    #K.P: This is a Archimedean spiral equation
    point = [i * np.cos(i), i * np.sin(i)]
    point_set.append(point)

max_radius = -1

for i in range(len(point_set)):
    point_minus_starting = [a - b for a,b in zip(point_set[i], viewpoint)]
    radius = alpha * np.linalg.norm(np.array(point_minus_starting))
    max_radius = max(radius, max_radius)

#K.P: Reflect the points across circle with radius R
reflected_points = []

for i in range(len(point_set)):
    cur_point = point_set[i]
    cur_norm = np.linalg.norm(np.array(cur_point))
    if cur_norm != 0:
        direction = cur_point / np.linalg.norm(np.array(cur_point))
    else:
        continue
    addition_vector = 2 * (max_radius - np.linalg.norm(cur_point))
    reflected_point = cur_point + direction * addition_vector
    reflected_dict[tuple(reflected_point)] = cur_point
    reflected_points.append(reflected_point)

#K.P: Add points to the plot 
plt.scatter(np.array(point_set)[:, 0], np.array(point_set)[:, 1], color='blue', s=invisible_thickness, label='Original Points')
plt.scatter(np.array(reflected_points)[:, 0], np.array(reflected_points)[:, 1], s=invisible_thickness, color='red', label='Reflected Points')

#K.P: Create the convex hull
convex_hull_input_points = reflected_points + point_set + [[0,0]] 
convex_hull = ConvexHull(convex_hull_input_points)

#K.P: Plot the Convex Hull 
for simplex in convex_hull.simplices:
    plt.plot(np.array(convex_hull_input_points)[simplex, 0], np.array(convex_hull_input_points)[simplex, 1], color="yellow")  

plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.title('Convex Hull Visualization')


#K.P: Definition 3: A point is visible if its inverted point lies on the convex hull
hull_vertices = np.array(convex_hull.points[convex_hull.vertices]) #K.P: Points making up the convex hull edges
visible_points = []

for point in reflected_points:
    #K.P: Point is visible if not a viewpoint and in the edge of the hull
    if np.any(np.all(point == hull_vertices, axis=1)) and not np.array_equal(point, viewpoint):
        visible_points.append(point)

#K.P: Darken the visible points and their corresponding actual values. 
for point in visible_points:
    plt.scatter(point[0], point[1], color='darkred', s=invisible_thickness * 2.2, edgecolors='red', label='Visible Points')  
    original_point = reflected_dict[tuple(point)]
    plt.scatter(original_point[0], original_point[1], color='darkblue', s=invisible_thickness * 2.2, edgecolors='blue', label='Visible Points')  

plt.show()