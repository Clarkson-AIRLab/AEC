% INPUT:1.U0 the local wind speed,(m/s),for aircraft U0=Uwind-Uaircraft
         %in the paper,the wind velocity in ms^-1 is the speed of
         %the surrounding air in relation to the sampling probe
       %2.U the flow velocity in the sampling prode(m/s)
       %3.Stk the Stokes Number of the sampling probe
       %4.theta_s the aspiration angle(o)
       %5.SOT the sampling orientiation, this parameter only has two value,
        %1=upword,0=downword
       %6.transmission define if transmission is considered or not, =1 considered =0 not considered
% OUTPUT: the transmission efficiency
function ef_trans=ef_trans_c(U0,U,Stk,theta_s,SOT,transmission)%calculate the trasmission efficiency
R=U0/U;
if transmission==1%consider the transmission efficiency
    if theta_s==0%isoaxial sampling
        if R<1%U0/U<1
            Iv=0.09.*(Stk.*(U-U0)/U0).^(0.3);
            ef_trans=exp(-75.*Iv.^2);%Hangal and Willeke (1990a,b)
%             if max(Stk)>4 || min(Stk)<0.02 || R<0.25%Validity range
%                  h = msgbox('waring:Stk must be in the range of [0.02,4] U0/U must in [0.5,2]', 'Error from transmission','error');
%            end
        elseif R>=1 && R<=10%1<=U0/U<=10
            ef_trans=(1+(U0/U-1)./(1+2.66.*Stk.^(-2/3)))./(1+(U0/U-1)./(1+0.418./Stk));%Liu et al.(1989)
%             if max(Stk)>100 || min(Stk)<0.01 || max(U0/U)>10 %Validity range
%                  h = msgbox('waring:Stk must be in the range of [0.01,100] U0/U must in [1,10]', 'Error from transmission','error');
%            end
        else
            ef_trans=1;
            warning('for transmission efficiency U0/U must be less than 10')
        end
    elseif theta_s>0 && theta_s <=90%non-isoaxial sampling 
        alpha=12*((1-theta_s/90)-exp(-theta_s));
        if SOT==1%upward
            Iw=Stk.*((U0/U)^(1/2))*sin((theta_s+alpha)*pi/180)*sin((theta_s+alpha)*pi/180/2);
        elseif SOT==0%downward
            Iw=Stk.*((U0/U)^(1/2))*sin((theta_s+alpha)*pi/180)*sin((theta_s-alpha)*pi/180/2);
        else
            warning('for the sampling orientiation, only the value 1 and 0 accept')
        end
        if R>=0.25 && R<=1
            Iv=0.09.*(Stk.*cos(theta_s*pi/180).*(U-U0)/U).^(0.3);
        else
            Iv=0;
        end
        ef_trans=exp((-75).*((Iv+Iw).^2));%Hangal and Willeke (1990a,b)
%         if max(Stk)>4 || min(Stk)<0.02 %Validity range
%                  h = msgbox('waring:Stk must be in the range of [0.02,4]', 'Error from transmission','error');
%         end
    else
        warning('for calculating the transmission efficiency, theta_s must be in the range of [0,90]')
    end
else
    ef_trans=1;
end
end