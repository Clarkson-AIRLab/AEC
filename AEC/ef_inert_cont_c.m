% INPUT:1.Re_tf the flow Re in the tubing 
       %2.Stk the Stokes Number of the sampling probe
       %3.theta_cont the contraction half-angle(o)
       %4.inertial bend define if bend is considered or not, =1 considered =0 not considered
       %5.Nt the number of the tubing section, Nt should always = length(L)
% OUTPUT: the efficiency of inertial contraction in tubing section
function ef_inert_cont=ef_inert_cont_c(dt,Stkt,theta_cont,contraction,Nt)%calculate the efficiency of inertial contraction
if contraction==1%consider the inertial deposition contraction efficiency Muyshondt et al.(1996b)
    counter=0;
   for i=1:Nt
       
%        if theta_cont(i)>90 || theta_cont(i) < 12%only consider the contracion for the range of the contration at [12,90]
%            ef_inert_cont(i,:)=1;
%        else
       Ai=dt(i)^2*pi;           %Ai is the cross-sectional area in front of the contraction,
       Ao=dt(i+1)^2*pi;         %Ao is the cross-sectional area behind the contraction.
       if Ai>Ao
           counter=counter+1;
           ef_inert_cont(counter,:)=1-1./(1+(Stkt(i,:).*(1-Ao/Ai)./(3.14*exp(-0.0185*theta_cont(i)))).^(-1.24));
%          cont_at(i)=i;
       end  
%        end
   end
%    if max(Stk)>100 || min(Stk)<0.001 %Validity range
%                  h = msgbox('waring:Stk must be in the range of [0.001,100]', 'Error from inertial deposition contraction','error');
%    end
else
    ef_inert_cont=1;
end