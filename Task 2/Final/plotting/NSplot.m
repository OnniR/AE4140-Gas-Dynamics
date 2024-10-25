function [] = NSplot(NSX,NSY, var, N)
    for j = 2:N
        for i = 1:N-1
            %if i+1 == j
            if i+1 == j && all([NSX(i,j), NSX(i+1,j), NSX(i+1,j+1)] ~= 0)
                avg = mean([var(i,j), var(i+1,j), var(i+1,j+1)]);
                plot_triangle(NSX(i,j), NSY(i,j), NSX(i+1,j), NSY(i+1,j), NSX(i+1,j+1), NSY(i+1,j+1), avg);
                plot_triangle(NSX(i,j), -NSY(i,j), NSX(i+1,j), -NSY(i+1,j), NSX(i+1,j+1), -NSY(i+1,j+1), avg)
            %else
            elseif all([NSX(i,j), NSX(i+1,j), NSX(i,j+1), NSX(i+1,j+1)] ~= 0)
                avg = mean([var(i,j), var(i+1,j), var(i, j+1), var(i+1,j+1)]);
                plot_square(NSX(i,j), NSY(i,j), NSX(i+1,j), NSY(i+1,j), NSX(i+1,j+1), NSY(i+1,j+1), NSX(i,j+1), NSY(i,j+1), avg)
                plot_square(NSX(i,j), -NSY(i,j), NSX(i+1,j), -NSY(i+1,j), NSX(i+1,j+1), -NSY(i+1,j+1), NSX(i,j+1), -NSY(i,j+1), avg)
            end
        end
    end
end