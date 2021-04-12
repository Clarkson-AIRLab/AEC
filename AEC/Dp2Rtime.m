%%%% function to calculate the relaxation time of a particle with diameter
%%%% Dp.
% INPUT: Dp in nm, Cc
% OUTPUT: relaxaion time

function tau = Dp2Rtime(Dp_nm,rho_p,K)
x=Dp_nm*1e-9; %changing to [m]
Cc=1 + (2.*K.mean_fp./x).*(K.alpha_Cc + K.beta_Cc.*exp(-K.gamma_Cc./(2.*K.mean_fp./x))); %Slip correction factor
tau=rho_p.*x.^2.*Cc/(18*K.mu); % relaxation time
end
