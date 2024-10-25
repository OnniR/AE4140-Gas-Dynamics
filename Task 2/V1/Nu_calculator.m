function Nu = Nu_calculator(nu1, phi1, phi2, plus_characteristic)
%Compute the Prandtl Meyer Angle along a characteristic line

if plus_characteristic
    Nu = phi2 - phi1 + nu1;
else 
    Nu = phi1 + nu1 - phi2;

end