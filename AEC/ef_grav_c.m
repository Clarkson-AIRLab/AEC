% INPUT:1.Re_tf the flow Re in the tubing 
       %2.Vts the terminal settling velocity of the particles,Willeke and Baron (2005)
       %3.Ut the average velocity in each tubing section(m/s)
       %4.L the length of each tubing section(m)
       %5.aver_dt the average inner diameter of tubing
       %6.theta_i the inclination(o)
       %7.sedimentation define if sedimentation is considered or not, =1 considered =0 not considered
       %8.Nt the number of the tubing section, Nt should always = length(L)
% OUTPUT: the efficiency of sedimentation in tubing section
function ef_grav=ef_grav_c(Re_tf,Vts,Ut,L,aver_dt,theta_i,sedimentation,Nt)%calculate the efficiency of sedimentation
if sedimentation==1%consider the diffusion efficiency
    for i=1:Nt%loops for different tubing sections
        if max(Vts)*sin(theta_i(i)*pi/180)/Ut(i) > 1
            warning('for calculating the sedimentation, the condition Vts*sin(theta_i)/Ut<<1 must satisify')
        end
        Z(i,:)=L(i).*Vts./(aver_dt(i)*Ut(i));%gravitational deposition parameter
        E(i,:)=3/4.*Z(i,:);
       if Re_tf(i)<2000%laminar flow
           if theta_i(i)==0%laminar flow in horizontal
               for j=1:length(Vts)
               if E(i,j)<=1
                   ef_grav(i,j)=1-2/pi*(2*E(i,j)*sqrt(1-E(i,j)^(2/3))-...
                   E(i,j)^(1/3)*sqrt(1-E(i,j)^(2/3))+asin(E(i,j)^(1/3)));%Fuchs 1964,Thomas 1958
               else
                   ef_grav(i,j)=0;
               end
               end
           else%laminar flow not in horizontal
               for j=1:length(Vts)
               Ek(i,j)=E(i,j).*cos(theta_i(i)*pi/180);
               if Ek(i,j)<=1
                   ef_grav(i,j)=1-2/pi*(2*Ek(i,j).*sqrt(1-Ek(i,j).^(2/3))-...
                   Ek(i,j).^(1/3).*sqrt(1-Ek(i,j).^(2/3))+asin(Ek(i,j).^(1/3)));% Heyder and Gebhurt 1977
               else
                   ef_grav(i,j)=0.*Ek(i,j);
               end
               end
           end
        else%turbulent flow Schwendiman et al. (1975)
           if theta_i(i)==0
               ef_grav(i,:)=exp(-4.*Z(i,:)./pi);
           else
               ef_grav(i,:)=exp(-4.*Z(i,:).*cos(theta_i(i)*pi/180)/pi);
           end 
       end
    end
else
    ef_grav=1;
end
end