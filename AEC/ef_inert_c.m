% INPUT:1.Re_tf the flow Re in the tubing 
       %2.Stk the Stokes Number of the sampling probe
       %3.Ut the average velocity in each tubing section(m/s)
       %4.L the length of each tubing section(m)
       %5.aver_dt the average inner diameter of tubing
       %6.Qt the flow rate in the tubing sections l/min
       %7.inertial define if inertial is considered or not, =1 considered =0 not considered
       %8.Nt the number of the tubing section, Nt should always = length(L)
% OUTPUT: the efficiency of inertial in tubing section
function ef_inert=ef_inert_c(Re_tf,Stkt,Ut,L,aver_dt,Qt,inertial,Nt)%calculate the efficiency of inertial
if inertial==1%consider the inertial deposition bend efficiency  Pui et al. 1987    
    for i=1:Nt%loops for different tubing sections
        Vt=((6*10^(-4).*(0.0395.*Stkt(i,:).*Re_tf(i).^(3/4)).^2+2*10^(-8).*Re_tf(i)).*Ut(i))./(5.03.*Re_tf(i).^(1/8));%willeke and Baron(2005)
        ef_inert(i,:)=exp(-pi*aver_dt(i)*L(i).*Vt/(Qt(i)*0.001/60));
    end
else
    ef_inert=1;
end
end