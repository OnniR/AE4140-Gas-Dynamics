function Nu = Nu_calculator(M, gamma)
% Nu_calculator: calculate the Prandtl Meyer Angle from a given Mach number

Nu = sqrt((gamma+1) / (gamma-1)) * atand(sqrt((gamma-1) / (gamma+1) * (M^2-1))) - atand(sqrt(M^2-1));

end