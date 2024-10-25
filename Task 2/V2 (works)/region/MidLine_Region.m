function [out] = MidLine_Region(S_Waves1T)

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
    NS_M(j-1, j) = M_calculator(NS_nu(j-1, j), gamma);
    MID_M(end+1) = NS_M(j-1, j);

    NS_mu(j-1, j) = Mu_calculator(NS_M(j-1, j));
    NS_alpha_minus(j-1, j) = NS_mu(j-1, j) - NS_phi(j-1, j);
    NS_alpha_plus(j-1, j) = NS_mu(j-1, j) + NS_phi(j-1, j);

    NS_x(j-1, j) = NS_y(j-1, j-1) / tand(NS_alpha_minus(j-1, j-1)) + NS_x(j-1, j-1);
    MID_x(end+1) = NS_x(j-1, j);

    NS_y(j-1, j) = 0;
    MID_y(end+1) = NS_y(j-1, j);
    
    if j < N+1
        for i = j:N
            NS_phi(i,j) = 0.5*(NS_phi(i,j-1) + NS_nu(i,j-1) + NS_phi(i-1,j) - NS_nu(i-1,j));
            NS_nu(i,j) = 0.5*(NS_phi(i,j-1) + NS_nu(i,j-1) + NS_nu(i-1,j) - NS_phi(i-1,j));
            NS_M(i,j) = M_calculator(NS_nu(i,j), gamma);
            NS_mu(i, j) = Mu_calculator(NS_M(i, j));
            NS_alpha_minus(i, j) = NS_mu(i, j) - NS_phi(i, j);
            NS_alpha_plus(i, j) = NS_mu(i, j) + NS_phi(i, j);
            NS_x(i,j) = (NS_y(i,j-1) + tand(NS_alpha_minus(i,j-1))*NS_x(i,j-1) - NS_y(i-1,j) + tand(NS_alpha_plus(i-1,j))*NS_x(i-1,j)) / (tand(NS_alpha_plus(i-1,j)) + tand(NS_alpha_minus(i,j-1)));
            NS_y(i,j) = tand(NS_alpha_plus(i-1,j))*NS_x(i,j) + NS_y(i-1,j) - tand(NS_alpha_plus(i-1,j))*NS_x(i-1,j);
        end
    end
out = [NS_phi(N,:)' NS_nu(N,:)' NS_M(N,:)' NS_mu(N,:)' NS_alpha_minus(N,:)' NS_alpha_plus(N,:)' NS_x(N,:)' NS_y(N,:)'];
end