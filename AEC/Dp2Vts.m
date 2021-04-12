%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: Diameter in nm and an array of constants (which includes the )
% OUTPUT: Settling velocity of the charged particle
function Vts = Dp2Vts(Dp_nm,rho_p,K)
x=Dp_nm*1e-9*sqrt(rho_p/1000); %changing to [m]
for i=1:length(x)
    Cc(i)=1 + (2.*K.mean_fp./x(i)).*(K.alpha_Cc + K.beta_Cc.*exp(-K.gamma_Cc/(2.*K.mean_fp./x(i)))); %Slip correction factor
    Vts(i)=rho_p.*(x(i)).^2.*K.g.*Cc(i)./(18*K.mu);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%