import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull


#K.P: Parameter initialization 
number_of_points: int = int(input("Enter number of points to generate: "))
alpha: int = 1.5
point_set= []
viewpoint_x= float(input("Enter a viewpoint X value:"))
viewpoint_y= float(input("Enter a viewpoint Y value:"))
viewpoint_z= -1 
dimensions: int = 2
reflected_dict = {}
invisible_thickness = 10
is_3D = (input("Are you using 3 dimensions? (Y/n)"))
is_3D = is_3D == "Y" or is_3D == "y"

def calc_max_radius(points, viewpoint):
    max_rad = -1
    for i in range(len(points)):
        point_minus_starting = [a - b for a,b in zip(points[i], viewpoint)]
        radius = alpha * np.linalg.norm(np.array(point_minus_starting))
        max_rad = max(radius, max_rad)
    return max_rad

if is_3D:
    viewpoint_z = float(input("Enter a viewpoint Z value:"))
    viewpoint = [viewpoint_x, viewpoint_y, viewpoint_z]
    dimensions = 3
    radius = 5  # Example radius
    center_x, center_y, center_z = 4, 0, 0  # Center of the sphere

    # Generate parameter values
    theta = np.linspace(0, 2 * np.pi, number_of_points)  # Azimuthal angle
    phi = np.linspace(0, np.pi, number_of_points)  # Polar angle

    point_set = []

    # Nested loops to cover the sphere surface
    for t in theta:
        for p in phi:
            x = radius * np.sin(p) * np.cos(t) + center_x
            y = radius * np.sin(p) * np.sin(t) + center_y
            z = radius * np.cos(p) + center_z
            point_set.append([x, y, z])
    
    max_radius = calc_max_radius(point_set, viewpoint)
    reflected_points = []
    for i in range(len(point_set)):
        cur_point = point_set[i]
        cur_norm = np.linalg.norm(np.array(cur_point))
        if cur_norm != 0:
            direction = cur_point / np.linalg.norm(np.array(cur_point))
        else:
            continue
        addition_vector = 2 * (max_radius - np.linalg.norm(cur_point))
        reflected_point = cur_point + (direction * addition_vector)
        reflected_dict[tuple(reflected_point)] = cur_point
        reflected_points.append(reflected_point)
    
    fig = plt.figure()  #
    ax = fig.add_subplot(111, projection='3d')  
    ax.scatter(np.array(point_set)[:, 0], np.array(point_set)[:, 1], np.array(point_set)[:, 2], color='blue', s=20, label='Original Points')

    # Reflected Points
    ax.scatter(np.array(reflected_points)[:, 0], np.array(reflected_points)[:, 1], np.array(reflected_points)[:, 2], color='red', s=20, label='Reflected Points')

    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    plt.legend()  

    convex_hull_input_points = reflected_points + point_set 
    convex_hull = ConvexHull(convex_hull_input_points)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    all_points = np.array(convex_hull_input_points)
    ax.scatter(all_points[:, 0], all_points[:, 1], all_points[:, 2], color="blue", s=20)

    # Plot the convex hull
    for simplex in convex_hull.simplices:
        # For 3D, simplices are indices of the vertices that form the faces of the convex hull
        # Therefore, we plot a line between each pair of vertices in each simplex (face)
        for i in range(len(simplex)):
            start_idx, end_idx = simplex[i], simplex[(i + 1) % len(simplex)]
            ax.plot(all_points[[start_idx, end_idx], 0], all_points[[start_idx, end_idx], 1], all_points[[start_idx, end_idx], 2], color="yellow")

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    plt.title('3D Convex Hull Visualization')
    
    # Prepare the 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Assuming hull_vertices is already defined as the vertices of the convex hull in 3D
    # And assuming reflected_dict maps 3D reflected points to their original 3D points

    visible_points = []
    hull_vertices = np.array(convex_hull.points[convex_hull.vertices])

    for point in reflected_points:
        # Check if the point is a vertex of the convex hull and not equal to the viewpoint
        if np.any(np.all(point == hull_vertices, axis=1)) and not np.array_equal(point, viewpoint):
            visible_points.append(point)

    # Plotting
    # Darken the visible points
    for point in visible_points:
        ax.scatter(point[0], point[1], point[2], color='darkred', s=40)  # Increase s for visibility
        original_point = reflected_dict[tuple(point)]
        ax.scatter(original_point[0], original_point[1], original_point[2], color='darkblue', s=40)  # Increase s for visibility

    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    plt.title('Visible Points in 3D')
    plt.show()
else:
    viewpoint = [viewpoint_x, viewpoint_y]
    #K.P: Generate point set
    xs = np.linspace(0, 5*np.pi, number_of_points)
    for i in xs:
        #K.P: This is the circle equation
        # point = [4 + np.cos(i), np.sin(i)]
        #K.P: This is a Archimedean spiral equation
        point = [i * np.cos(i) + 20  , i * np.sin(i) + 0]
        point_set.append(point)

    max_radius = calc_max_radius(point_set, viewpoint)

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
        reflected_point = cur_point + (direction * addition_vector)
        reflected_dict[tuple(reflected_point)] = cur_point
        reflected_points.append(reflected_point)

    #K.P: Add points to the plot 
    plt.scatter(np.array(point_set)[:, 0], np.array(point_set)[:, 1], color='blue', s=invisible_thickness, label='Original Points')
    plt.scatter(np.array(reflected_points)[:, 0], np.array(reflected_points)[:, 1], s=invisible_thickness, color='red', label='Reflected Points')

    #K.P: Create the convex hull
    convex_hull_input_points = reflected_points + point_set 
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
        plt.scatter(point[0], point[1], color='darkred', s=invisible_thickness * 2.2)  
        original_point = reflected_dict[tuple(point)]
        plt.scatter(original_point[0], original_point[1], color='darkblue', s=invisible_thickness * 2.2)  

    plt.show()



