% INPUT:1.the cunningham slip correction factor(for oil droplets and solid particles),
         %Allen and Raabe,(1985)
       %2.the aerodynamic particle diameter(nm), Willeke and Baron,(2005)
       %3.Ut the flow velocity inside the tubing section
       %4.%the particle density (kg/m^3)
       %5.aver_dt the average diameter in each tubing section
       %6.% the air dynamic viscosity (Pa*s) or(kg/(m*s))
% OUTPUT: the Stokes Number of the each tubing sections
function Stkt=Stkt_c(Cc,Dp_aero,Ut,rho_p,aver_dt,K)
da=Dp_aero.*10^(-9);
for i=1:length(aver_dt)
Stkt(i,:)=da.^2.*rho_p.*Cc.*Ut(i)./(18.*K.mu.*aver_dt(i));
end
end