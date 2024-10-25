
function plot_square(x1, y1, x2, y2, x3, y3, x4, y4, M)
    % X and Y coordinates of the square
    xv = [x1, x2, x3, x4];
    yv = [y1, y2, y3, y4];

    % Add edge line for better visibility
    hold on;
    %plot([xv, xv(1)], [yv, yv(1)], 'k-');  % To close the square

    % Fill the square with a constant color based on the Mach number M
    fill(xv, yv, M, "LineStyle","none");

    % Axis labels
    xlabel('X');
    ylabel('Y');
end