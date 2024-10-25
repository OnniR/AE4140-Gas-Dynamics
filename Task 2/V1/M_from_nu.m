function M = M_from_nu(nu, gamma)
% Nu_calculator: calculate the Prandtl Meyer Angle from a given Mach number

syms M
nu_eq = Nu_from_M(M, gamma) == nu;
M = abs(vpasolve(nu_eq, M));

end