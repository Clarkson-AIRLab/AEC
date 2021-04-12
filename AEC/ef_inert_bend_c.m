% INPUT:1.Re_tf the flow Re in the tubing 
       %2.Stk the Stokes Number of the sampling probe
       %3.theta_Kr the angle of curvature(o)
       %4.inertial bend define if bend is considered or not, =1 considered =0 not considered
       %5.Nt the number of the tubing section, Nt should always = length(L)
% OUTPUT: the efficiency of inertial bend in tubing section
function ef_inert_bend=ef_inert_bend_c(Re_tf,Stkt,theta_Kr,bend,Nt)%calculate the efficiency of inertial bend
if bend==1%consider the inertial deposition bend efficiency  Pui et al. 1987
    for i=1:Nt%loops for different tubing sections
        if Re_tf(i)<2000%laminar flow
                if theta_Kr(i)==0
                   ef_inert_bend(i,:)=exp(Stkt(i,:).*theta_Kr(i)*pi/180); 
                else
                   ef_inert_bend(i,:)=(1+(Stkt(i,:)./0.171).^(0.452.*(Stkt(i,:)./0.171+2.24))).^(-2/pi/(theta_Kr(i)*pi/180));
                end
            else%turbulent flow
                ef_inert_bend(i,:)=exp(-2.823.*Stkt(i,:).*theta_Kr(i)*pi/180); 
        end
    end
else
    ef_inert_bend=1;
end
end