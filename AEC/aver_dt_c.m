% INPUT:1.dt the determined diameter of the tubing,(m) in the paper,
        %the first diameter of tubing should be the same as the d0
       %2.Nt the number of the tubing section, Nt should always = length(L)
% OUTPUT: aver_dt the average inner diameter of tubing
function aver_dt=aver_dt_c(area_aver,Nt)
for i=1:Nt
    aver_dt(i)=2*sqrt(area_aver(i)/pi);%the average inner diameter of tubing
end
end