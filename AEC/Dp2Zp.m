%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: Diameter in nm and an array of constants (which includes the )
% OUTPUT: Mobility of the charged particle
function [Zp] = Dp2Zp(Dp_nm,charge,K)
x=Dp_nm*1e-9; %changing to [m]
for i=1:length(x)
    Cc(i)=1 + (2*K.mean_fp./x(i)).*(K.alpha_Cc + K.beta_Cc*exp(-K.gamma_Cc./(2*K.mean_fp./x(i)))) ;%Slip correction factor
    Zp(i)=(charge.*K.e.*Cc(i)./(3.*pi.*K.mu*(x(i))));  %m^2.V^-1.s^-1 
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%