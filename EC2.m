%K.P: Make the meshgrid of the unit square [0,1]^2 with intervals of 0.25
[X, Y] = meshgrid(0:0.05:1, 0:0.05:1);

xv = X(:);
yv = Y(:);

%K.P: Use Delaunay to make a triangulation from the X Y meshgrid
DT = delaunayTriangulation(xv,yv);


%K.P: Visualize the triangulation with triplot
triplot(DT);

%K.P: Find boundary points
boundaryPoints = freeBoundary(DT);

%K.P: Initialize the boundary condition array
z = zeros(size(xv));

bpts = [0,1];

%K.P: Apply sin(pi * y) for all X boundary points and -sin(pi*x) for all
%Y values
for i = 1:size(boundaryPoints, 1)
    pointIndex = boundaryPoints(i, 1);
    x = xv(pointIndex);
    y = yv(pointIndex);
    if ismember(x, bpts)
        z(pointIndex) = sin(pi*y);  
    end
    
    pointIndex = boundaryPoints(i, 2);
    y = yv(pointIndex);
    if ismember(y, bpts)
        z(pointIndex) = -sin(pi*x);  
    end
end
%% 

cvx_begin
    variable Z(length(xv))
    %K.P: Apply boundary conditions
    Z(boundaryPoints) == z(boundaryPoints);
    
    %K.P: Minimize the areas of the triangles
    %K.P: connectivityList will give the triangle indicies

    expression areas(size(DT.ConnectivityList, 1), 1);
    for i = 1:size(DT.ConnectivityList, 1)
        verts = DT.ConnectivityList(i, :);

        %K.P: Point 1
        X1 = xv(verts(1)); 
        Y1 = yv(verts(1)); 
        Z1 = Z(verts(1));
        %K.P: Point 2
        X2 = xv(verts(2)); 
        Y2 = yv(verts(2)); 
        Z2 = Z(verts(2));
        %K.P: Point 3
        X3 = xv(verts(3)); 
        Y3 = yv(verts(3)); 
        Z3 = Z(verts(3));
        
        %K.P: Vector from vertex 1 to vertex 2
        V12 = [X2-X1; Y2-Y1; Z2-Z1];
        % Vector from vertex 1 to vertex 3
        V13 = [X3-X1; Y3-Y1; Z3-Z1];
        
        %K.P: Cross product of V12 and V13
        cross_prod = cross(V12, V13);
        
        %K.P: Calculate the area of the triangle 
        areas(i) = norm(cross_prod, 2) / 2;
        
    end
    minimize(sum(areas))
cvx_end

%% 

%K.P: Get the triangle verticies
tri = DT.ConnectivityList;

%K.P: Plot the shape
trimesh(tri, xv, yv, Z);  % 'k' sets the edge color to black
axis equal;  
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Optimized Mesh Shape with Visible Triangles');
shading interp;  
colorbar;  
view(3); 