function [] = Splot(NSX,NSY,var,N)
for i = 1:N-1
    j = 1;
    avg = mean([var(i,j), var(i,j+1), var(i+1,j), var(i+1,j+1)]);
    plot_square(NSX(i,j), NSY(i,j), NSX(i+1,j), NSY(i+1,j), NSX(i+1,j+1),  NSY(i+1,j+1), NSX(i,j+1), NSY(i,j+1), avg);
    plot_square(NSX(i,j), -NSY(i,j), NSX(i+1,j), -NSY(i+1,j), NSX(i+1,j+1),  -NSY(i+1,j+1), NSX(i,j+1), -NSY(i,j+1), avg);
end
end