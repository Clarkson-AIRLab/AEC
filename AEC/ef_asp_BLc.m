% INPUT:1.U0 the local wind speed,(m/s),for aircraft U0=Uwind-Uaircraft
         %in the paper,the wind velocity in ms^-1 is the speed of
         %the surrounding air in relation to the sampling probe
       %2.U the flow velocity in the sampling prode(m/s)
       %3.Stk the Stokes Number of the sampling probe
       %4.theta_s the aspiration angle(o)
       %5.theta_v the angle corresponding to the vertical(o)
       %6.rho_p the particle density
       %7.the aerodynamic particle diameter(nm)
       %8.aspiration define if aspiration is considered or not, =1 considered =0 not considered
% OUTPUT: the aspiration efficiency
function ef_asp=ef_asp_BLc(U0,U,Stk,theta_s,theta_v,rho_p,Dp_nm,Dp_aero,K,Cc,aspiration)
da=Dp_aero.*10^(-9);
X=1;%the shape factor
g=K.g;
mu=K.mu;
dp=Dp_nm*10^-9;
if aspiration==1%consider the aspiration efficiency
    if U0>1.5%for moving air
        %isoaxial sampling condition
        if theta_s==0
            if  U0/U>5.6 && U0/U<=10
                k=2+06.2*(U/U0)-0.91*(U/U0)^0.103;%Paik and Vincent (2002)
            elseif U0/U<=5.6
                k=2+0.317*(U/U0);
            else
                warning('For B&L sampler U0/U cannot higher than 10')
            end
           ef_asp=1+(U0/U-1)*(1-1./(1+k.*Stk));% Belyaev and Levin (1972,1974)
%            if max(Stk)>2.03 || min(Stk)<0.05 ||U0/U<0.17%Validity range, Stevens (1986)
%                  h = msgbox('waring:Stk must be in the range of [0.05,2.03]', 'Error from Aspiration','error');
%            end
        %Non-isoaxial sampling condition
        elseif theta_s>0 && theta_s<=60
            Stkk=Stk.*exp(0.022*theta_s);
            ef_asp=1+(U0/U*cos(theta_s*pi/180)-1)*...
                ((1-(1+(2+0.617*(U/U0)).*Stkk).^(-1))/(1-(1+2.617.*Stkk)))*...
                (1-(1 + 0.55.*Stkk.*exp(0.25.*Stkk)).^(-1));%Durham and Lundgren (1980)
%             if max(Stk)>4 || min(Stk)<0.02 || U0/U>2 || U0/U<0.5%Validity range
%                  h = msgbox('waring:Stk must be in the range of [0.02,4] U0/U must in [0.5,2]', 'Error from Aspiration','error');
%            end
        elseif theta_s>=61 && theta_s<=90
            ef_asp=1+(U0/U*cos(theta_s*pi/180)-1)*(3.*Stk.^(sqrt(U0/U)));%Hangal and Willeke (1990a)
%             if max(Stk)>0.2 || min(Stk)<0.02 || U0/U>2 || U0/U<0.5%Validity range
%                  h = msgbox('waring:Stk must be in the range of [0.02,0.2] U0/U must in [0.5,2]', 'Error from Aspiration','error');
%             end
        else
            warning('For moving air input theta_s must be in the range of [0,90]')%this will stop the program
        end
    elseif U0<0.5%calm air
             Vts=(rho_p.*da.^2*g.*Cc)/(18*X*mu);%Willeke and Baron (2005)
            ef_asp=Vts./U*cos(theta_v*pi/180)+exp(-(4*(Stk.^(1+(Vts./U).^(1/2)))/(1+2.*Stk)));%Grinshpun et al. (1993 1994)
           Rep=Vts.*rho_p.*dp./mu;%the particle Re
%            if max(Stk)>100 || min(Stk)<0.001 || max(Vts./U)>2 || min(Vts./U)<0.5%Validity range
%                  h = msgbox('waring:Stk must be in the range of [0.001,100] U0/U must in [0.5,2]', 'Error from Aspiration','error');
%            end
           if  theta_v<0 ||theta_v>90
               warning('for calm air input theta_v must be in the range of [0,90]')
           end
           if max(Rep)>=0.1
               warning('for calm air the particle Re must be less than 0.1');
           end
    else%slow motion 0.5<=U0<=1.5
        V0=U0;%V0 is the initial velocity of the particle,here set the V0=the wind velocity around the samping probe
        Vts=V0-U0;
        phi=(Vts./U0)*(Vts./U0+2*cos((theta_s+theta_v)*pi/180));
        f_moving=exp(-Vts./U0);
        f_calm=1-exp(-Vts./U0);
         %calculate the ef_asp_moving
         if theta_s==0
           k=2+0.317*(U/U0);
           ef_asp_moving=1+(U0/U-1)*(1-1./(1+k.*Stk));% Belyaev and Levin (1972,1974)
        elseif theta_s>0 || theta_s<=60
            Stkk=Stk.*exp(0.022*theta_s);
            ef_asp_moving=1+(U0./U*cos(theta_s*pi/180)-1)*...
                ((1-(1+(2+0.617*(U./U0)).*Stkk).^(-1))/(1-(1+2.617.*Stkk)))*...
                (1-(1 + 0.55.*Stkk.*exp(0.25.*Stkk)).^(-1));%Durham and Lundgren (1980)
        elseif theta_s>=61 && theta_s<=90
            ef_asp_moving=1+(U0./U*cos(theta_s*pi/180)-1)*(3*Stk^(sqrt(U0./U)));%Hangal and Willeke (1990a)
        else
            warning('for moving air input theta_s must be in the range of [0,90]')
        end
         %calculate the ef_asp_calm
        Vts_calm=(rho_p.*da.^2*g.*Cc)/(18*X*mu);%Willeke and Baron (2005)
        ef_asp_calm=Vts_calm./U*cos(theta_v)+exp(-(4*(Stk.^(1+sqrt(Vts_calm./U)))/(1+2.*Stk)));%Grinshpun et al. (1993 1994)
        if  theta_v<0 ||theta_v>90
            warning('for calm air input theta_v must be in the range of [0,90], and the particle Re should <0.1')
        end
        %ef_asp for slow motion
        ef_asp=ef_asp_moving.*sqrt(1+phi)*f_moving+ef_asp_calm.*f_calm;%Grinshpun et al. (1993, 1994)
    end
else
    ef_asp=1;
end
end