% INPUT:1.Re_tf the flow Re in the tubing 
       %2.D the particle diffusion coefficient D
       %3.L the length of each tubing section(m)
       %4.Qt the flow rate in the tubing sections l/min
       %5.Sc the Schmidt Number Sc
       %6.diffusion define if diffusion is considered or not, =1 considered =0 not considered
       %7.Nt the number of the tubing section, Nt should always = length(L)
% OUTPUT: the efficiency of diffusion in tubing section
function ef_diffusion=ef_diffusion_c(Re_tf,D,L,Qt,Sc,diffusion,Nt)%calculate the efficiency of diffusion
if diffusion==1%consider the diffusion efficiency
    for i=1:Nt%loops for different tubing sections
        Xi(i,:)=pi.*D.*L(i)./(Qt(i)*0.001/60);
        if Re_tf(i)<2000%laminar flow
            Sh(i,:)=3.66+0.2672./(Xi(i,:)+0.10079.*Xi(i,:).^(0.3));%the Sherwood Number, Holman (1972)
        else%turbulent flow
            Sh(i,:)=0.0118*Re_tf(i).^(7/8)*Sc.^(1/3);%the Sherwood Number, Friedlander and Johnstone(1957)
        end
       ef_diffusion(i,:)=exp(-Xi(i,:).*Sh(i,:));%Willeke and Buron (2005)
    end
else
    ef_diffusion=1;    
end
end