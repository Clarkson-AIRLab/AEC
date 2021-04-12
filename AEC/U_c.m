% INPUT:1.d0_mm the orifice diameter (mm),in the paper 
         %the orifice diameter defind as the inner diameter of the tube opening, 
         %at the point where the aerosol enters the tubing,
         %so we use d0 to calculate the flow velocity in samping prode
       %2.Q the flow rate in sampling prode(L/min)
         %in the paper,The flow rate in l min?1 is 
         %that measured in the first tube section immediately downstream of the orifice
% OUTPUT: the flow velocity in the sampling prode
function U=U_c(d0_mm,Q)
d0 = d0_mm*1e-3;
U= Q*0.001/(60*pi*(d0/2)^2);%(m/s)
end