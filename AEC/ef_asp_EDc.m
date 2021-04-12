% INPUT:1.U0 the local wind speed,(m/s),for aircraft U0=Uwind-Uaircraft
         %in the paper,the wind velocity in ms^-1 is the speed of
         %the surrounding air in relation to the sampling probe
       %2.U the flow velocity in the sampling prode(m/s)
       %3.Stk the Stokes Number of the sampling probe
       %4.rho_f the flow density of inside sampler
       %5.rho_0 the flow density of outside air
       %6.d0_mm the inner diameter(mm)of the sampler
       %7.dout_mm the outter diameter(mm) of the sampler
       %8.dbb_m, a spherical shaped blund body diameter
       %9.L_inlet_m, a distance from the inlet entrance
       %10.aspiration define if aspiration is considered or not, =1 considered =0 not considered
% OUTPUT: the aspiration efficiency
function ef_asp=ef_asp_EDc(U0,U,Stk,rho_0,rho_f,d0_mm,dout_mm,dbb_mm,L_inlet_mm,aspiration)%calculate the aspiration efficien
din=d0_mm*10^(-3);
dout=dout_mm*10^(-3);
dbb=dbb_mm*10^(-3);
Lin=L_inlet_mm*10^(-3);
if aspiration==1%consider the aspiration efficiency
    R3=rho_0*U0/(rho_f*U);
    R1=0.02*R3+1.19;
    R2=0.78*(dout/din)+0.17*(dbb/Lin)+0.52*(dout/din)*(dbb/Lin)-0.59;
    ef_asp=(R3-1)./(1+R2.*Stk.^(-R1))+1;
else
    ef_asp=1;
end
end