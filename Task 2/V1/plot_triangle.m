function plot_triangle(x1, y1, x2, y2, x3, y3, M)
    % X and Y coordinates of the triangle
    xv = [x1, x2, x3];
    yv = [y1, y2, y3];
    
    
    % Add edge line for better visibility
    hold on;
    %plot([xv, xv(1)], [yv, yv(1)], 'k-');  % To close the triangle

    % Fill the triangle with a constant color based on the Mach number M
    fill(xv, yv, M, "LineStyle","none");
    
    % Axis labels
    %xlabel('X');
    %ylabel('Y');
end
