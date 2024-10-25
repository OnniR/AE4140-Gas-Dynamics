
function plot_square(x1, y1, x2, y2, x3, y3, x4, y4, M)
    % X and Y coordinates of the square
    xv = [x1, x2, x3, x4];
    yv = [y1, y2, y3, y4];

    hold on;

    % Fill the square with a constant color based on the Mach number M
    %fill(xv, yv, M, "LineStyle","none", 'HandleVisibility', 'off'); %Plot without characteristic
    fill(xv, yv, "w", 'HandleVisibility', 'off'); %Plot with characteristic lines in black and no fill

end