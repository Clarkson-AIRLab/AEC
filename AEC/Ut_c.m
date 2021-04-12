% INPUT:1.Qt%the flow rate in the tubing sections l/min
       %2.Nt the number of the tubing section, Nt should always = length(L)
% OUTPUT: Ut the average velocity in each tubing section(m/s)
function Ut=Ut_c(Qt,area_aver,Nt)
for i=1:Nt
    Ut(i)=Qt(i)*0.001/60/area_aver(i);%average velocity in each tubing section(m/s)
end
end