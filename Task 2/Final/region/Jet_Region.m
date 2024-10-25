function [NS2_x, NS2_y, NS2_phi, NS2_nu, NS2_M, NS2_mu, NS2_alpha_minus, NS2_alpha_plus, JB_x, JB_y, JB_phi] = Jet_Region(out, N, gamma, JB_x, JB_y, JB_phi, nu_2)
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
NS2_phi(:,1) = out(:,1);
NS2_nu(:,1) = out(:,2);
NS2_M(:,1) = out(:,3);
NS2_mu(:,1) = out(:,4);
NS2_alpha_minus(:,1) = out(:,5);
NS2_alpha_plus(:,1) = out(:,6);
NS2_x(:,1) = out(:,7);
NS2_y(:,1) = out(:,8);


for j = 2:N+1
    % Compute values on the jet boundary
    NS2_nu(j-1, j) = nu_2;
    NS2_phi(j-1, j) = NS2_phi(j-1,j-1) - NS2_nu(j-1,j-1) + NS2_nu(j-1,j);

    NS2_M(j-1, j) = M_calculator(NS2_nu(j-1, j), gamma);
    NS2_mu(j-1, j) = Mu_calculator(NS2_M(j-1, j));
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
            NS2_M(i,j) = M_calculator(NS2_nu(i,j), gamma);
            NS2_mu(i, j) = Mu_calculator(NS2_M(i, j));
            NS2_alpha_minus(i, j) = NS2_mu(i, j) - NS2_phi(i, j);
            NS2_alpha_plus(i, j) = NS2_mu(i, j) + NS2_phi(i, j);

            NS2_x(i,j) = (-NS2_y(i,j-1) + tand(NS2_alpha_plus(i,j-1))*NS2_x(i,j-1) + NS2_y(i-1,j) + tand(NS2_alpha_minus(i-1,j))*NS2_x(i-1,j)) / (tand(NS2_alpha_minus(i-1,j)) + tand(NS2_alpha_plus(i,j-1)));
            NS2_y(i,j) = tand(NS2_alpha_plus(i,j-1))*NS2_x(i,j) - tand(NS2_alpha_plus(i,j-1))*NS2_x(i,j-1) + NS2_y(i,j-1);

        end
    end

end
end