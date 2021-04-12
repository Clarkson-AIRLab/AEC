%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: Zp and charge
% OUTPUT: Dp_nm
function [Dp_star1] = Zp2Dp(Zp_star,charge,K)%diameter[m],charge,pressure,temp
% x=1e-9*(1:1000);
%P in kPa
if length(size(Zp_star))==2
[row,col]=size(Zp_star);depth=1;
elseif length(size(Zp_star))==3
   [row,col,depth]=size(Zp_star);
end

for i_r=1:row
    for i_c=1:col
        for i_d=1:depth
    range_Zp=1e-12;
    dDp=1;
    Cc=1;
    while(abs(dDp)>range_Zp)
        Dp_star=(charge.*K.e*Cc./(3.*pi.*K.mu*Zp_star(i_r,i_c,i_d)));
        Cc= 1 + (2.*K.mean_fp/Dp_star).*(K.alpha_Cc + K.beta_Cc.*exp(-K.gamma_Cc...
            ./(2.*K.mean_fp./Dp_star))); %Slip correction factor
        Dp_star2=(charge.*K.e.*Cc./(3.*pi.*K.mu*Zp_star(i_r,i_c,i_d)));
        dDp=(Dp_star2-Dp_star)./Dp_star2;
        Dp_star=Dp_star2;
    end
    Dp_star1(i_r,i_c,i_d)=Dp_star;
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%