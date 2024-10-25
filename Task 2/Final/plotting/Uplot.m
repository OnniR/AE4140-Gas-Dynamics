function [] = Uplot(NSX, NSY, NS2X, NS2Y, var, var2, N)
avg = mean([var(N,1) var2(1,1), var2(1,2)]);
plot_triangle(NSX(N,1), NSY(N,1), NS2X(1,1), NS2Y(1,1), NS2X(1,2), NS2Y(1,2), avg);
plot_triangle(NSX(N,1), -NSY(N,1), NS2X(1,1), -NS2Y(1,1), NS2X(1,2), -NS2Y(1,2), avg);

end
