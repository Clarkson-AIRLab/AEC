% INPUT: Dp_nm
% OUTPUT: aerodynamic diameter ignoring the slip correction

function Dp_aero = Dp2Dp_aero(Dp_nm,rho_p)
Dp_aero=Dp_nm*sqrt(rho_p/1000);

end