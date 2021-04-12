%% INPUT: Dp_arr
%% OUTPUT: CPC efficiency

function [cpc_eff]=cpc_efficiency(Dp_nm,type)

switch type
    case 'TSI_3010'
        M=[6.4 9.9];
    case 'TSI_3025'
        M=[1.9 3.6];
    case 'TSI_3786'
        M=[3.2 3.9];
end

Dp=1:0.01:1000;
eta=zeros(1,length(Dp));
% TSI_3010=[6.4 9.9];
% TSI_3025=[1.9 3.6];
% TSI_3786=[3.2 3.9];
% M=[TSI_3010;
% TSI_3025;
% TSI_3786];
% M=TSI_3010;
% for i=1:size(M,1)
    D0=M(1,1);
    D50=M(1,2);
    alp1=D0;
    alp2=(D50-D0)/log(2);
    for j=1:length(Dp)
        if Dp(j)>=D0
            eta(:,j)=1-exp((alp1-Dp(j))./alp2);
        else
            eta(:,j)=0;
        end
    end
%     semilogx(Dp_arr,eta)
%     hold on
% end
% hold off
cpc_eff=interp1(Dp,eta,Dp_nm);
end