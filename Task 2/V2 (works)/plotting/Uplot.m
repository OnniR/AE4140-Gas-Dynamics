function [] = Uplot(Boundary_x, Boundary_y, Region_x, Region_y, var, N, region_index)
boundary_index = 1 + (region_index-1)*N
plot_triangle(Boundary_x(boundary_index),Boundary_y(boundary_index), Region_x(region_index,region_index),Region_y(region_index,region_index),Boundary_x(boundary_index+1),Boundary_y(boundary_index+1), var(region_index,region_index));
plot_triangle(Boundary_x(boundary_index),-Boundary_y(boundary_index), Region_x(region_index,region_index),-Region_y(region_index,region_index),Boundary_x(boundary_index+1),-Boundary_y(boundary_index+1), var(region_index,region_index));

end
