clear all;
close all;
clc;

%________________Initial Values________________%
gamma = 1.4;
pa = 100000;
pe = 2.5*pa;
Me = 2;
phi_e = 0;
H = 1;
N = 10;

%________________Initial Computations________________%
%Compute the flow properties of the uniform region ACD with the initial conditions in the nozzle
M2 = M_IsentropicRelation(pe, pa, Me, gamma); %Calculate M_ACD from the isentropic relations
nu_e = Nu_calculator(Me, gamma); %Compute nu_e 
nu_2 = Nu_calculator(M2, gamma); %Compute nu_ACD
plus = true;
phi2 = phi_e - nu_e + nu_2; %Compute phi in region ACD using a Gamma+ characteristic from region e

%________________Jet Boundary Data________________%
%Initialise vectors to store the flow direction and coordinates of reflections points of characteristics on the jet boundary
JB_phi = [phi2]; 
JB_x = [0];
JB_y = [H];

%________________Mid-line Data________________%
%Initialise vectors to store the coordinates of reflections points of characteristics on the center line
MID_x = [0];
MID_y = [0];

%________________Wave Divider________________%
D_phi = (phi2 - phi_e) / (N - 1); % Divide the total expansion wave into N equally distributed waves
phi_arr = [];j
for i = 1:N % Compute the angle of each wave
    phi_i = D_phi*(i-1);
    phi_arr = [phi_arr, phi_i];
end

%________________First Simple Expansion Fan________________%
%Initialise vectors to store flow properties on each characteristic line
nu_arr = [nu_e];
M_arr = [];
mu_arr = [];
x_arr = [];
y_arr = [];

%Given the intial value for and nu and all the known values for phi from the wave divider, we can compute nu on the next expansion wave along a Gamma+ characteristic
for i = 2:N 
    nu_i = phi_arr(i) - phi_arr(i-1) + nu_arr(i-1); 
    nu_arr = [nu_arr, nu_i];
end

%Compute M, mu, the coordinates of point A, alpha+ and alpha- and store them in a vector
for i = 1:length(nu_arr);
    M_i = M_calculator(nu_arr(i), gamma);
    M_arr = [M_arr, M_i];
    mu_i = Mu_calculator(M_arr(i));
    mu_arr = [mu_arr, double(mu_i)];
    x_arr = [x_arr, 0];
    y_arr = [y_arr, H];
end
alpha_min_arr = mu_arr - phi_arr;
alpha_plus_arr = mu_arr + phi_arr;

%Combine all data vectors in a single matrix and round off for a cleaner look when printing the matrix 
init_data = round([phi_arr; nu_arr; M_arr; mu_arr; alpha_min_arr; alpha_plus_arr; x_arr; y_arr], 5); %Columns correspond to charactersitics, rows to flow properties

%________________Following Regions________________%

[NS_x, NS_y, NS_phi, NS_nu, NS_M, NS_mu, NS_alpha_minus, NS_alpha_plus, MID_x, MID_y] = MidLine_Region(init_data', N, gamma, MID_x, MID_y);
out = [NS_phi(N,2:N+1)' NS_nu(N,2:N+1)' NS_M(N,2:N+1)' NS_mu(N,2:N+1)' NS_alpha_minus(N,2:N+1)' NS_alpha_plus(N,2:N+1)' NS_x(N,2:N+1)' NS_y(N,2:N+1)'];

[NS2_x, NS2_y, NS2_phi, NS2_nu, NS2_M, NS2_mu, NS2_alpha_minus, NS2_alpha_plus, JB_x, JB_y, JB_phi] = Jet_Region(out, N, gamma, JB_x, JB_y, JB_phi, nu_2);
out2 = [NS2_phi(N,2:N+1)' NS2_nu(N,2:N+1)' NS2_M(N,2:N+1)' NS2_mu(N,2:N+1)' NS2_alpha_minus(N,2:N+1)' NS2_alpha_plus(N,2:N+1)' NS2_x(N,2:N+1)' NS2_y(N,2:N+1)'];

[NS3_x, NS3_y, NS3_phi, NS3_nu, NS3_M, NS3_mu, NS3_alpha_minus, NS3_alpha_plus, MID_x, MID_y] = MidLine_Region(out2, N, gamma, MID_x, MID_y);
out3 = [NS3_phi(N,2:N+1)' NS3_nu(N,2:N+1)' NS3_M(N,2:N+1)' NS3_mu(N,2:N+1)' NS3_alpha_minus(N,2:N+1)' NS3_alpha_plus(N,2:N+1)' NS3_x(N,2:N+1)' NS3_y(N,2:N+1)'];


%________________Plotting________________%
hold on

%Plot jet boundary
plot(JB_x, JB_y, "--", "Color","black")
plot(JB_x, -JB_y, "--", "Color","black")
legend("Jet Boundary")

%Plot jet exhaust region
Splot(NS_x, NS_y, NS_M, N) %Plot first simple region
NSplot(NS_x, NS_y, NS_M, N) %Plot first non-simple region
plot_triangle(0,0,0,H,NS_x(1,2),NS_y(1,2),Me) %Plot first uniform region
plot_triangle(0,0,0,-H,NS_x(1,2),-NS_y(1,2),Me) %Plot first uniform region
Splot(NS2_x, NS2_y, NS2_M, N) %Plot second simple region
NSplot(NS2_x, NS2_y, NS2_M, N) %Plot second non-simple region
Uplot(NS_x,NS_y,NS2_x,NS2_y,NS_M,NS2_M,N) %Plot second uniform region
Splot(NS3_x, NS3_y, NS3_M, N) %Plot second simple region
NSplot(NS3_x, NS3_y, NS3_M, N) %Plot second non-simple region
Uplot(NS2_x,NS2_y,NS3_x,NS3_y,NS2_M,NS3_M,N) %Plot second uniform region


axis equal
colormap("parula")
colorbarHandle = colorbar()
caxis([Me max(max(NS2_M))])
xlabel("X")
ylabel("Y")
ylabel(colorbarHandle, 'Mach Number')
grid on
grid minor










%% ________________First Non-Simple Wave________________%
% 
% %Initiate NxN+1 matrix for each variable of interest
% NS_phi = zeros(N, N+1);
% NS_nu = zeros(N, N+1);
% NS_M = zeros(N, N+1);
% NS_mu = zeros(N, N+1);
% NS_alpha_minus = zeros(N, N+1);
% NS_alpha_plus = zeros(N, N+1);
% NS_x = zeros(N, N+1);
% NS_y = zeros(N, N+1);
% 
% %Replace first column with data from simple region
% NS_phi(:,1) = S_Waves1T(:,1);
% NS_nu(:,1) = S_Waves1T(:,2);
% NS_M(:,1) = S_Waves1T(:,3);
% NS_mu(:,1) = S_Waves1T(:,4);
% NS_alpha_minus(:,1) = S_Waves1T(:,5);
% NS_alpha_plus(:,1) = S_Waves1T(:,6);
% NS_x(:,1) = 0;
% NS_y(:,1) = H;
% 
% 
% for j = 2:N+1
%     % Compute values on the center-line of the jet
%     NS_phi(j-1, j) = 0;
%     NS_nu(j-1, j) = NS_phi(j-1, j-1) + NS_nu(j-1,j-1);
%     NS_M(j-1, j) = M_calculator(NS_nu(j-1, j), gamma);
%     MID_M(end+1) = NS_M(j-1, j);
% 
%     NS_mu(j-1, j) = Mu_calculator(NS_M(j-1, j));
%     NS_alpha_minus(j-1, j) = NS_mu(j-1, j) - NS_phi(j-1, j);
%     NS_alpha_plus(j-1, j) = NS_mu(j-1, j) + NS_phi(j-1, j);
% 
%     NS_x(j-1, j) = NS_y(j-1, j-1) / tand(NS_alpha_minus(j-1, j-1)) + NS_x(j-1, j-1);
%     MID_x(end+1) = NS_x(j-1, j);
% 
%     NS_y(j-1, j) = 0;
%     MID_y(end+1) = NS_y(j-1, j);
% 
%     if j < N+1
%         for i = j:N
%             NS_phi(i,j) = 0.5*(NS_phi(i,j-1) + NS_nu(i,j-1) + NS_phi(i-1,j) - NS_nu(i-1,j));
%             NS_nu(i,j) = 0.5*(NS_phi(i,j-1) + NS_nu(i,j-1) + NS_nu(i-1,j) - NS_phi(i-1,j));
%             NS_M(i,j) = M_calculator(NS_nu(i,j), gamma);
%             NS_mu(i, j) = Mu_calculator(NS_M(i, j));
%             NS_alpha_minus(i, j) = NS_mu(i, j) - NS_phi(i, j);
%             NS_alpha_plus(i, j) = NS_mu(i, j) + NS_phi(i, j);
%             NS_x(i,j) = (NS_y(i,j-1) + tand(NS_alpha_minus(i,j-1))*NS_x(i,j-1) - NS_y(i-1,j) + tand(NS_alpha_plus(i-1,j))*NS_x(i-1,j)) / (tand(NS_alpha_plus(i-1,j)) + tand(NS_alpha_minus(i,j-1)));
%             NS_y(i,j) = tand(NS_alpha_plus(i-1,j))*NS_x(i,j) + NS_y(i-1,j) - tand(NS_alpha_plus(i-1,j))*NS_x(i-1,j);
%         end
%     end
% 
% end
% 
% 
%% %________________Second Non-Simple Wave________________%
% 
% %Initiate NxN+1 matrix for each variable of interest
% NS2_phi = zeros(N, N+1);
% NS2_nu = zeros(N, N+1);
% NS2_M = zeros(N, N+1);
% NS2_mu = zeros(N, N+1);
% NS2_alpha_minus = zeros(N, N+1);
% NS2_alpha_plus = zeros(N, N+1);
% NS2_x = zeros(N, N+1);
% NS2_y = zeros(N, N+1);
% 
% %Replace first column with data from simple region
% NS2_phi(:,1) = NS_phi(N,2:end);
% NS2_nu(:,1) = NS_nu(N,2:end);
% NS2_M(:,1) = NS_M(N,2:end);
% NS2_mu(:,1) = NS_mu(N,2:end);
% NS2_alpha_minus(:,1) = NS_alpha_minus(N,2:end);
% NS2_alpha_plus(:,1) = NS_alpha_plus(N,2:end);
% NS2_x(:,1) = NS_x(N,2:end);
% NS2_y(:,1) = NS_y(N,2:end);
% 
% 
% for j = 2:N+1
%     % Compute values on the jet boundary
%     NS2_nu(j-1, j) = nu_2;
%     NS2_phi(j-1, j) = NS2_phi(j-1,j-1) - NS2_nu(j-1,j-1) + NS2_nu(j-1,j);
% 
%     NS2_M(j-1, j) = M_calculator(NS_nu(j-1, j), gamma);
%     JB_M(end+1) = NS_M(j-1, j);
% 
%     NS2_mu(j-1, j) = Mu_calculator(NS_M(j-1, j));
%     NS2_alpha_minus(j-1, j) = NS2_mu(j-1, j) - NS2_phi(j-1, j);
%     NS2_alpha_plus(j-1, j) = NS2_mu(j-1, j) + NS2_phi(j-1, j);
% 
%     % Compute coordinates of the jet boundary points
%     JB_x(end+1) = (JB_y(end) - JB_x(end)*tand(JB_phi(end)) - NS2_y(j-1,j-1) + NS2_x(j-1,j-1)*tand(NS2_alpha_plus(j-1,j-1))) / (-tand(JB_phi(end)) + tand(NS2_alpha_plus(j-1,j-1)));
%     NS2_x(j-1, j) = JB_x(end);
% 
%     JB_y(end+1) = tand(JB_phi(end))*JB_x(end) + JB_y(end) - tand(JB_phi(end))*JB_x(end-1);
%     NS2_y(j-1, j) = JB_y(end);
% 
%     % Append jet boundary phi's to array
%     JB_phi(end+1) = NS2_phi(j-1,j);
% 
% 
%     if j < N+1
%         for i = j:N
%             NS2_phi(i,j) = 0.5*(NS2_phi(i-1,j) + NS2_nu(i-1,j) + NS2_phi(i,j-1) - NS2_nu(i,j-1));
%             NS2_nu(i,j) = 0.5*(NS2_phi(i-1,j) + NS2_nu(i-1,j) - NS2_phi(i,j-1) + NS2_nu(i,j-1));
%             NS2_M(i,j) = M_calculator(NS2_nu(i,j), gamma);
%             NS2_mu(i, j) = Mu_calculator(NS_M(i, j));
%             NS2_alpha_minus(i, j) = NS2_mu(i, j) - NS2_phi(i, j);
%             NS2_alpha_plus(i, j) = NS2_mu(i, j) + NS2_phi(i, j);
% 
%             NS2_x(i,j) = (-NS2_y(i,j-1) + tand(NS2_alpha_plus(i,j-1))*NS2_x(i,j-1) + NS2_y(i-1,j) + tand(NS2_alpha_minus(i-1,j))*NS2_x(i-1,j)) / (tand(NS2_alpha_minus(i-1,j)) + tand(NS2_alpha_plus(i,j-1)));
%             NS2_y(i,j) = tand(NS2_alpha_plus(i,j-1))*NS2_x(i,j) - tand(NS2_alpha_plus(i,j-1))*NS2_x(i,j-1) + NS2_y(i,j-1);
% 
%         end
%     end
% 
% end




%% Plot first simple region
% for i = 1:N-1
%     j = 2;
%     avg = mean([NS_M(i,j), NS_M(i+1,j), NS_M(i+1,j-1)]);
%     plot_triangle(NS_x(i,j), NS_y(i,j), NS_x(i+1,j), NS_y(i+1,j), NS_x(i+1,j-1), NS_y(i+1,j-1), avg)
% end
%% Plot first non simple region
% for j = 2:N
%     for i = 1:N-1
%         if i+1 == j
%             avg = mean([NS_M(i,j), NS_M(i+1,j), NS_M(i+1,j+1)]);
%             plot_triangle(NS_x(i,j), NS_y(i,j), NS_x(i+1,j), NS_y(i+1,j), NS_x(i+1,j+1), NS_y(i+1,j+1), avg)
%         else
%             avg = mean([NS_M(i,j), NS_M(i+1,j), NS_M(i, j+1), NS_M(i+1,j+1)]);
%             plot_square(NS_x(i,j), NS_y(i,j), NS_x(i+1,j), NS_y(i+1,j), NS_x(i+1,j+1), NS_y(i+1,j+1), NS_x(i,j+1), NS_y(i,j+1), avg)
%         end
%     end
% end
%% Plot second simple region
% for i = 1:N-1
%     j = 1;
%     avg = mean([NS2_M(i,j), NS2_M(i,j+1), NS2_M(i+1,j), NS2_M(i+1,j+1)]);
%     plot_square(NS2_x(i,j), NS2_y(i,j), NS2_x(i+1,j), NS2_y(i+1,j), NS2_x(i+1,j+1),  NS2_y(i+1,j+1), NS2_x(i,j+1), NS2_y(i,j+1), avg);
% end
%% Plot second non simple region
% for j = 2:N
%     for i = 1:N-1
%         if i+1 == j
%             avg = mean([NS2_M(i,j), NS2_M(i+1,j), NS2_M(i+1,j+1)]);
%             plot_triangle(NS2_x(i,j), NS2_y(i,j), NS2_x(i+1,j), NS2_y(i+1,j), NS2_x(i+1,j+1), NS2_y(i+1,j+1), avg)
%         else
%             avg = mean([NS2_M(i,j), NS2_M(i+1,j), NS2_M(i, j+1), NS2_M(i+1,j+1)]);
%             plot_square(NS2_x(i,j), NS2_y(i,j), NS2_x(i+1,j), NS2_y(i+1,j), NS2_x(i+1,j+1), NS2_y(i+1,j+1), NS2_x(i,j+1), NS2_y(i,j+1), avg)
%         end
%     end
% end

