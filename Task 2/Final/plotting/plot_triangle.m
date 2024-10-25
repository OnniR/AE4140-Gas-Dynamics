function plot_triangle(x1, y1, x2, y2, x3, y3, M)
    % X and Y coordinates of the triangle
    xv = [x1, x2, x3];
    yv = [y1, y2, y3];
    
    hold on;

    % Fill the triangle with a constant color based on the Mach number M
    %fill(xv, yv, M, "LineStyle","none", 'HandleVisibility', 'off'); %Plot without characteristic lines
    fill(xv, yv, "w", 'HandleVisibility', 'off'); %Plot with characteristic lines in black and no fill
end
