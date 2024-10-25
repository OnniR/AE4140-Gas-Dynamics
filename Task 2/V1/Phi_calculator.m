function Phi = Phi_calculator(nu1, nu2, phi1, plus_characteristic)
%Phi_calculator: calculates phi using the invariant of either a plus or a minus characteristic

if plus_characteristic
    Phi = phi1 - nu1 + nu2;
else 
    Phi = phi1 + nu1 - nu2;

end