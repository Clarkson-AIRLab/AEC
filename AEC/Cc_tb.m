% INPUT: diameter in nm,K
% OUTPUT: Slip correction
function [Cc]=Cc_tb(Dp,K)
Dp=Dp*1e-9;
Cc= 1 + (2.*K.mean_fp./Dp).*(K.alpha_Cc + K.beta_Cc.*exp(-K.gamma_Cc ...
            ./(2.*K.mean_fp./Dp))); %Slip correction factor
end

