function M2 = M_calculator(p1, p2, M1, gamma)
%M_calculator: compute the mach number from the isentropic relations

M2 = sqrt(2/(gamma-1) * ((1 + (gamma-1)/2*M1^2) * (p1/p2)^((gamma-1)/gamma)-1));

end