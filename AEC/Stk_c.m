% INPUT:1.the cunningham slip correction factor(for oil droplets and solid particles),
         %Allen and Raabe,(1985)
       %2.the aerodynamic particle diameter(nm), Willeke and Baron,(2005)
       %3.U0 the local wind speed,(m/s),for aircraft U0=Uwind-Uaircraft
         %in the paper,the wind velocity in ms^-1 is the speed of
         %the surrounding air in relation to the sampling probe
       %4.%the particle density (kg/m^3)
       %5.%the orifice diameter (mm),in the paper 
         %the orifice diameter defind as the inner diameter of the tube opening, 
         %at the point where the aerosol enters the tubing,
         %so we use d0 to calculate the flow velocity in samping prode
       %6.% the air dynamic viscosity (Pa*s) or(kg/(m*s))
% OUTPUT: the Stokes Number of the sampling probe
function Stk=Stk_c(Cc,Dp_aero,U0,rho_p,d0_mm,K)
d0 = d0_mm*1e-3;
da=Dp_aero.*10^(-9);
Stk=da.^2.*rho_p.*Cc.*U0./(18.*K.mu.*d0);
end