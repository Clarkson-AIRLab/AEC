% INPUT:1.dt the determined diameter of the tubing,(m) in the paper,
        %the first diameter of tubing should be the same as the d0
       %2.Nt the number of the tubing section, Nt should always = length(L)
% OUTPUT: area_aver the average area for each tubing section 
function area_aver=area_aver_c(dt,Nt)
for i=1:Nt
    area_aver(i)=(pi/3)*((dt(i+1)/2)^2+dt(i)*dt(i+1)/4+(dt(i)/2)^2);%the average area of truncated cone
end
end