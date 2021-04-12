function [ef_dma,mu_dma]=ef_dma_c(D,Qa,effL_cm)
effL=effL_cm*10^(-2);
for i=1:length(D)

    mu_dma(i)=D(i)*effL/(Qa*0.001/60);

if mu_dma(i)< 0.009
    ef_dma(i)=1-5.5*mu_dma(i)^0.666 + 3.77*mu_dma(i);
else
    ef_dma(i)=0.819*exp(-11.5*mu_dma(i)) + 0.0975*exp(-70.1*mu_dma(i));
end
end

