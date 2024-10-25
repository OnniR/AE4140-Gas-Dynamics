clear all;
close all;
clc;

%________________Initial Values________________%
N = 3;
gamma = 1.4;
pa = 100000;
pe = 2*pa;
Me = 2;
phi_e = 0;
H = 1;
N = 4;

%________________Initial Computations________________%
M2 = M_calculator(pe, pa, Me, gamma);
nu_e = Nu_from_M(Me, gamma);
nu_2 = Nu_from_M(M2, gamma);
plus = true;
phi2 = Phi_calculator(nu_e, nu_2, phi_e, plus);

%________________Jet Boundary Data________________%
JB_phi = [phi2]; %jet boundary phi's
JB_x = [0];
JB_y = [H];


%________________Wave Divider________________%
D_phi = (phi2 - phi_e) / (N - 1); % Divide the total expansion wave into N equally distributed waves
phi_arr = [];
for i = 1:N % Compute the angle of each wave
    phi_i = D_phi*(i-1);
    phi_arr = [phi_arr, phi_i];
end


%________________Simple Wave________________%
nu_arr = [nu_e];
plus = true;
for i = 2:N
    nu_i = Nu_calculator(nu_arr(i-1), phi_arr(i-1), phi_arr(i), plus);
    nu_arr = [nu_arr, nu_i];
end

M_arr = [];
for i = 1:length(nu_arr);
    M_i = M_from_nu(nu_arr(i), gamma);
    M_arr = [M_arr, M_i];
end

mu_arr = [];
for i = 1:length(M_arr);
    mu_i = Mu_calculator(M_arr(i));
    mu_arr = [mu_arr, double(mu_i)];
end

alpha_min_arr = mu_arr - phi_arr;
alpha_plus_arr = mu_arr + phi_arr;

S_Waves1 = round([phi_arr; nu_arr; M_arr; mu_arr; alpha_min_arr; alpha_plus_arr], 5); %Columns correspond to charactersitics, rows to flow properties
S_Waves1T = transpose(S_Waves1); 


%________________First Non-Simple Wave________________%

%Initiate NxN+1 matrix for each variable of interest
NS_phi = zeros(N, N+1);
NS_nu = zeros(N, N+1);
NS_M = zeros(N, N+1);
NS_mu = zeros(N, N+1);
NS_alpha_minus = zeros(N, N+1);
NS_alpha_plus = zeros(N, N+1);
NS_x = zeros(N, N+1);
NS_y = zeros(N, N+1);

%Replace first column with data from simple region
NS_phi(:,1) = S_Waves1T(:,1);
NS_nu(:,1) = S_Waves1T(:,2);
NS_M(:,1) = S_Waves1T(:,3);
NS_mu(:,1) = S_Waves1T(:,4);
NS_alpha_minus(:,1) = S_Waves1T(:,5);
NS_alpha_plus(:,1) = S_Waves1T(:,6);
NS_x(:,1) = 0;
NS_y(:,1) = H;


for j = 2:N+1
    % Compute values on the center-line of the jet
    NS_phi(j-1, j) = 0;
    NS_nu(j-1, j) = NS_phi(j-1, j-1) + NS_nu(j-1,j-1);
    NS_M(j-1, j) = M_from_nu(NS_nu(j-1, j), gamma);
    NS_mu(j-1, j) = Mu_calculator(NS_M(j-1, j));
    NS_alpha_minus(j-1, j) = NS_mu(j-1, j) - NS_phi(j-1, j);
    NS_alpha_plus(j-1, j) = NS_mu(j-1, j) + NS_phi(j-1, j);
    NS_x(j-1, j) = NS_y(j-1, j-1) / tand(NS_alpha_minus(j-1, j-1)) + NS_x(j-1, j-1);
    NS_y(j-1, j) = 0;

    
    if j < N+1
        for i = j:N
            NS_phi(i,j) = 0.5*(NS_phi(i,j-1) + NS_nu(i,j-1) + NS_phi(i-1,j) - NS_nu(i-1,j));
            NS_nu(i,j) = 0.5*(NS_phi(i,j-1) + NS_nu(i,j-1) + NS_nu(i-1,j) - NS_phi(i-1,j));
            NS_M(i,j) = M_from_nu(NS_nu(i,j), gamma);
            NS_mu(i, j) = Mu_calculator(NS_M(i, j));
            NS_alpha_minus(i, j) = NS_mu(i, j) - NS_phi(i, j);
            NS_alpha_plus(i, j) = NS_mu(i, j) + NS_phi(i, j);
            NS_x(i,j) = (NS_y(i,j-1) + tand(NS_alpha_minus(i,j-1))*NS_x(i,j-1) - NS_y(i-1,j) + tand(NS_alpha_plus(i-1,j))*NS_x(i-1,j)) / (tand(NS_alpha_plus(i-1,j)) + tand(NS_alpha_minus(i,j-1)));
            NS_y(i,j) = tand(NS_alpha_plus(i-1,j))*NS_x(i,j) + NS_y(i-1,j) - tand(NS_alpha_plus(i-1,j))*NS_x(i-1,j);
        end
    end

end


%________________Second Non-Simple Wave________________%

%Initiate NxN+1 matrix for each variable of interest
NS2_phi = zeros(N, N+1);
NS2_nu = zeros(N, N+1);
NS2_M = zeros(N, N+1);
NS2_mu = zeros(N, N+1);
NS2_alpha_minus = zeros(N, N+1);
NS2_alpha_plus = zeros(N, N+1);
NS2_x = zeros(N, N+1);
NS2_y = zeros(N, N+1);

%Replace first column with data from simple region
NS2_phi(:,1) = NS_phi(N,2:end);
NS2_nu(:,1) = NS_nu(N,2:end);
NS2_M(:,1) = NS_M(N,2:end);
NS2_mu(:,1) = NS_mu(N,2:end);
NS2_alpha_minus(:,1) = NS_alpha_minus(N,2:end);
NS2_alpha_plus(:,1) = NS_alpha_plus(N,2:end);
NS2_x(:,1) = NS_x(N,2:end);
NS2_y(:,1) = NS_y(N,2:end);


for j = 2:N+1
    % Compute values on the jet boundary
    NS2_nu(j-1, j) = nu_2;
    NS2_phi(j-1, j) = NS2_phi(j-1,j-1) - NS2_nu(j-1,j-1) + NS2_nu(j-1,j);

    NS2_M(j-1, j) = M_from_nu(NS_nu(j-1, j), gamma);
    NS2_mu(j-1, j) = Mu_calculator(NS_M(j-1, j));
    NS2_alpha_minus(j-1, j) = NS2_mu(j-1, j) - NS2_phi(j-1, j);
    NS2_alpha_plus(j-1, j) = NS2_mu(j-1, j) + NS2_phi(j-1, j);

    % Compute coordinates of the jet boundary points
    JB_x(end+1) = (JB_y(end) - JB_x(end)*tand(JB_phi(end)) - NS2_y(j-1,j-1) + NS2_x(j-1,j-1)*tand(NS2_alpha_plus(j-1,j-1))) / (-tand(JB_phi(end)) + tand(NS2_alpha_plus(j-1,j-1)));
    NS2_x(j-1, j) = JB_x(end);

    JB_y(end+1) = tand(JB_phi(end))*JB_x(end) + JB_y(end) - tand(JB_phi(end))*JB_x(end-1);
    NS2_y(j-1, j) = JB_y(end);
        
    % Append jet boundary phi's to array
    JB_phi(end+1) = NS2_phi(j-1,j);


    if j < N+1
        for i = j:N
            NS2_phi(i,j) = 0.5*(NS2_phi(i-1,j) + NS2_nu(i-1,j) + NS2_phi(i,j-1) - NS2_nu(i,j-1));
            NS2_nu(i,j) = 0.5*(NS2_phi(i-1,j) + NS2_nu(i-1,j) - NS2_phi(i,j-1) + NS2_nu(i,j-1));
            NS2_M(i,j) = M_from_nu(NS2_nu(i,j), gamma);
            NS2_mu(i, j) = Mu_calculator(NS_M(i, j));
            NS2_alpha_minus(i, j) = NS2_mu(i, j) - NS2_phi(i, j);
            NS2_alpha_plus(i, j) = NS2_mu(i, j) + NS2_phi(i, j);

            NS2_x(i,j) = (-NS2_y(i,j-1) + tand(NS2_alpha_plus(i,j-1))*NS2_x(i,j-1) + NS2_y(i-1,j) + tand(NS2_alpha_minus(i-1,j))*NS2_x(i-1,j)) / (tand(NS2_alpha_minus(i-1,j)) + tand(NS2_alpha_plus(i,j-1)));
            NS2_y(i,j) = tand(NS2_alpha_plus(i,j-1))*NS2_x(i,j) - tand(NS2_alpha_plus(i,j-1))*NS2_x(i,j-1) + NS2_y(i,j-1);

        end
    end

end


%________________Plotting________________%
hold on
plot(NS_x, NS_y, "o")
plot(NS2_x, NS2_y, "o")

%plot(NS2_x, NS2_y, "o")

plot(NS_x, NS_y, 'o')

%Plot jet exhaust region
plot_triangle(0,0, 0,H, NS_x(1,2),NS_y(1,2), 2)
%Plot first simple region
for i = 1:N-1
    j = 2;
    avg = mean([NS_M(i,j), NS_M(i+1,j), NS_M(i+1,j-1)]);
    plot_triangle(NS_x(i,j), NS_y(i,j), NS_x(i+1,j), NS_y(i+1,j), NS_x(i+1,j-1), NS_y(i+1,j-1), avg)
end
%Plot first non simple region
for j = 2:N
    for i = 1:N-1
        if i+1 == j
            avg = mean([NS_M(i,j), NS_M(i+1,j), NS_M(i+1,j+1)]);
            plot_triangle(NS_x(i,j), NS_y(i,j), NS_x(i+1,j), NS_y(i+1,j), NS_x(i+1,j+1), NS_y(i+1,j+1), avg)
        else
            avg = mean([NS_M(i,j), NS_M(i+1,j), NS_M(i, j+1), NS_M(i+1,j+1)]);
            plot_square(NS_x(i,j), NS_y(i,j), NS_x(i+1,j), NS_y(i+1,j), NS_x(i+1,j+1), NS_y(i+1,j+1), NS_x(i,j+1), NS_y(i,j+1), avg)
        end
    end
end
%Plot second non simple region
for j = 2:N
    for i = 1:N-1
        if i+1 == j
            avg = mean([NS2_M(i,j), NS2_M(i+1,j), NS2_M(i+1,j+1)]);
            plot_triangle(NS2_x(i,j), NS2_y(i,j), NS2_x(i+1,j), NS2_y(i+1,j), NS2_x(i+1,j+1), NS2_y(i+1,j+1), avg)
        else
            avg = mean([NS2_M(i,j), NS2_M(i+1,j), NS2_M(i, j+1), NS2_M(i+1,j+1)]);
            plot_square(NS2_x(i,j), NS2_y(i,j), NS2_x(i+1,j), NS2_y(i+1,j), NS2_x(i+1,j+1), NS2_y(i+1,j+1), NS2_x(i,j+1), NS2_y(i,j+1), avg)
        end
    end
end

axis equal
colormap("parula")
colorbar()
caxis([Me max(max(NS_M))])
xlim([0 max(max(NS2_x))])
ylim([0 max(max(NS2_y))])

NS2_x;



