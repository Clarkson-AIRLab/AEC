function [Properties,Input,Others]=AEC(varargin)

clc;
% the config.txt file contains a list of all the default inputs like
% temperature, pressure, charge, density. It also contains a list of
% constants commonly used in the aerosol toolbox like the Boltzsman const.,
% alpha, beta, gamma to calculate the slip correction, acceleration due to
% gravity, 1 ev in Joules etc. 
% the config file is loaded to access all the default inputs and constants
% to perform operations
fid=fopen('config.txt');
scan=textscan(fid,'%s %f','delimiter','=');
fclose(fid);

% Depending on what the user inputs, the default values are updated. The
% varargin is an array of texts, ex: varargin={'Dp_nm',100,'charge',1,'Zp',1e-8}
% using 'find' to locate where a particular variable is inputted, the
% corresponding value assigned to the variable is updated. 

variables={'charge','T','rho_p','P','alpha_Cc','beta_Cc','gamma_Cc','air_molecule_diameter','rho_0','rho_in'};
for i_check=1:length(variables)
ind(i_check)=isempty(find(strcmpi(varargin,variables{i_check})==1,1));
end

ind_update=find(ind==0);

if(nargin ==1 ||  nargin ==0  )
    disp('Error: Not enough inputs')
    return
elseif(isempty(find(strcmpi(varargin,'Dp_nm')==1, 1))*isempty(find(strcmpi(varargin,'Zp')==1, 1))==1)
    disp('Specify either Dp_nm or Zp')
    return
elseif(isempty(find(strcmpi(varargin,'Dp_nm')==1, 1))==0 && isempty(find(strcmpi(varargin,'Zp')==1, 1))==0)
    disp('Error: can''t specify Dp_nm and Zp')
    return
elseif(isempty(find(strcmpi(varargin,'Dp_nm')==1, 1))==0 && isempty(find(strcmpi(varargin,'Zp')==1, 1))==1)
    Dp_nm=varargin{find(strcmpi(varargin,'Dp_nm')==1 )+1};
elseif(isempty(find(strcmpi(varargin,'Dp_nm')==1, 1))==1 && isempty(find(strcmpi(varargin,'Zp')==1, 1))==0)
    Zp=varargin{find(strcmpi(varargin,'Zp')==1 )+1};
end


for i_up=1:length(variables)
inputs(i_up)=scan{2}(strcmpi(scan{1},variables{i_up}));
end

%%%% Update the inputs based on varargin
for i_up=1:length(ind_update)
ind_var(i_up)=  find(strcmpi(varargin,variables{ind_update(i_up)})==1); % find the position of the variable that's to be updated from varargin
inputs(ind_update(i_up))=varargin{ind_var(i_up)+1};
end

charge=inputs(strcmpi(variables,'charge'));
K.T=inputs(strcmpi(variables,'T'));
rho_p=inputs(strcmpi(variables,'rho_p'));
rho_0=inputs(strcmpi(variables,'rho_0'));
rho_f=inputs(strcmpi(variables,'rho_in'));
K.P=inputs(strcmpi(variables,'P'));
K.alpha_Cc=inputs(strcmpi(variables,'alpha_Cc'));
K.beta_Cc=inputs(strcmpi(variables,'beta_Cc'));
K.gamma_Cc=inputs(strcmpi(variables,'gamma_Cc'));
K.air_molecule_diameter=inputs(strcmpi(variables,'air_molecule_diameter'));
K.k=scan{2}(strcmpi(scan{1},'k'));
K.e=scan{2}(strcmpi(scan{1},'e'));
K.g=scan{2}(strcmpi(scan{1},'g'));
K.mean_fp=K.k*K.T/(sqrt(2)*pi*K.air_molecule_diameter^2*K.P); 
K.mu = 1.51204*K.T^1.5/(K.T+120)*1e-6   ;              % air dynamic viscosity in Pa s , changes for other gases like nitrogen/ oxygen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(isempty(find(strcmpi(varargin,'Zp'), 1))==1)
[Zp] = Dp2Zp(Dp_nm,charge,K);
elseif(isempty(find(strcmpi(varargin,'Dp_nm')==1, 1)))
    Dp = Zp2Dp(Zp,charge,K);
    Dp_nm=Dp*1e9;
end

Cc=Cc_tb(Dp_nm,K);
Vts = Dp2Vts(Dp_nm,rho_p,K);
tau_p = Dp2Rtime(Dp_nm,rho_p,K);
% S=
Dp_aero = Dp2Dp_aero(Dp_nm,rho_p);%add the standard density, here use water density=1000 (Kg/(m^3))
Kn=2*K.mean_fp*1e9./Dp_nm;
%calculate the particle diffusion coefficient D
D=K.k*K.T.*Cc./(3*pi*K.mu.*Dp_nm*10^(-9));
%calculate the Schmidt Number Sc
Sc=K.mu/rho_f./D;

%%%%%%  aerosol code use for effeciency analysis   %%%%%
%effeciency of sampling probe
%the input variables
efs=fopen('sampling.txt');
sampling_scan=textscan(efs,'%s %f','delimiter','=');
fclose(efs);
sampling_variables={'U0','Q','d0_mm','SOT','theta_s','theta_v','dout_mm','dbb_mm','L_inlet_mm',};
for i_check_in=1:length(sampling_variables)
ind_in(i_check_in)=isempty(find(strcmpi(varargin,sampling_variables{i_check_in})==1,1));
end

ind_update_in=find(ind_in==0);

for i_up=1:length(sampling_variables)
sampling_inputs(i_up)=sampling_scan{2}(strcmpi(sampling_scan{1},sampling_variables{i_up}));
end

for i_up=1:length(ind_update_in)
ind_var_in(i_up)=  find(strcmpi(varargin,sampling_variables{ind_update_in(i_up)})==1);
sampling_inputs(ind_update_in(i_up))=varargin{ind_var_in(i_up)+1};
end

U0=sampling_inputs(strcmpi(sampling_variables,'U0'));
Q=sampling_inputs(strcmpi(sampling_variables,'Q'));
d0_mm=sampling_inputs(strcmpi(sampling_variables,'d0_mm'));
SOT=sampling_inputs(strcmpi(sampling_variables,'SOT'));
theta_s=sampling_inputs(strcmpi(sampling_variables,'theta_s'));
theta_v=sampling_inputs(strcmpi(sampling_variables,'theta_v'));
dout_mm=sampling_inputs(strcmpi(sampling_variables,'dout_mm'));
dbb_mm=sampling_inputs(strcmpi(sampling_variables,'dbb_mm'));
L_inlet_mm=sampling_inputs(strcmpi(sampling_variables,'L_inlet_mm'));

%calculation
U=U_c(d0_mm,Q);%the flow velocity in sampling probe
Stk=Stk_c(Cc,Dp_aero,U0,rho_p,d0_mm,K);%the Stokes Number of the sampling probe
Re=rho_f*U*d0_mm*0.001/K.mu;%the flow Re in the sampling probe
if(5-(isempty(find(strcmpi(varargin,'B&L')==1, 1))+ isempty(find(strcmpi(varargin,'Eddy')==1, 1))...
        +isempty(find(strcmpi(varargin,'SMAI5')==1, 1))+isempty(find(strcmpi(varargin,'SMAI20')==1, 1))...
        +isempty(find(strcmpi(varargin,'BASE')==1, 1))) >=2)
    disp('Error: can not specify sampler type, type in one sampler each time')
    return
elseif isempty(find(strcmpi(varargin,'B&L')==1, 1))==0 
    aspiration=1;transmission=1;
    ef_asp=ef_asp_BLc(U0,U,Stk,theta_s,theta_v,rho_p,Dp_nm,Dp_aero,K,Cc,aspiration);%calculate the aspiration efficiency
    ef_trans=ef_trans_c(U0,U,Stk,theta_s,SOT,transmission);%calculate the trasmission efficiency
    %calculate the overall sampling efficiency
    ef_sampling=ef_asp.*ef_trans;%Willeke and Baron (2005)
    %semilogx(Dp_nm,ef_sampling)%for test
elseif isempty(find(strcmpi(varargin,'Eddy')==1, 1))==0
    aspiration=1;transmission=1;
    ef_asp=ef_asp_EDc(U0,U,Stk,rho_0,rho_f,d0_mm,dout_mm,dbb_mm,L_inlet_mm,aspiration);%calculate the aspiration efficiency
    ef_trans=ef_trans_c(U0,U,Stk,theta_s,SOT,transmission);%calculate the trasmission efficiency
    %calculate the overall sampling efficiency
    ef_sampling=ef_asp.*ef_trans;%Willeke and Baron (2005)
    %semilogx(Stk,ef_asp)%for test
elseif isempty(find(strcmpi(varargin,'SMAI5')==1, 1))==0
    ef_asp=1;
    ef_trans=1;
    [dp,ef] = importfile('SMAI5.txt');
    for i=1:length(Dp_nm)
    if Dp_nm(i)*10^(-3)>=0.01 && Dp_nm(i)*10^(-3)<=10
          ef_sampling(i)=interp1(dp,ef,Dp_nm(i)*10^(-3));
    elseif Dp_nm(i)*10^(-3)<0.01
        ef_sampling(i)=1;
    elseif Dp_nm(i)*10^(-3)>10
        ef_sampling(i)=0;
    end
    end
    %semilogx(Dp_nm,ef_sampling)%for test
elseif isempty(find(strcmpi(varargin,'SMAI20')==1, 1))==0
    ef_asp=1;
    ef_trans=1;
    [dp,ef] = importfile('SMAI20.txt');
    for i=1:length(Dp_nm)
    if Dp_nm(i)*10^(-3)>=0.01 && Dp_nm(i)*10^(-3)<=4
          ef_sampling(i)=interp1(dp,ef,Dp_nm(i)*10^(-3));
    elseif Dp_nm(i)*10^(-3)<0.01
        ef_sampling(i)=1;
    elseif Dp_nm(i)*10^(-3)>4
        ef_sampling(i)=0;
    end
    end
    %semilogx(Dp_nm,ef_sampling)%for test
elseif isempty(find(strcmpi(varargin,'BASE')==1, 1))==0
    ef_asp=1;
    ef_trans=1;
    [dp,ef] = importfile('BASEstreamline.txt');
    for i=1:length(Dp_nm)
    if Dp_nm(i)*10^(-3)>=0.01 && Dp_nm(i)*10^(-3)<=3
          ef_sampling(i)=interp1(dp,ef,Dp_nm(i)*10^(-3));
    elseif Dp_nm(i)*10^(-3)<0.01
        ef_sampling(i)=1;
    elseif Dp_nm(i)*10^(-3)>3
        ef_sampling(i)=0;
    end
    end
    %semilogx(Dp_nm,ef_sampling)%for test
else
    ef_asp=1;ef_trans=1;ef_sampling=1;
end


%effeciency of partical loss in tubing
%the input variables
if isempty(find(strcmpi(varargin,'dt_mm'), 1))==0
    dt_mm=varargin{find(strcmpi(varargin,'dt_mm')==1 )+1};
    if dt_mm(1)==d0_mm
        dt_mm=varargin{find(strcmpi(varargin,'dt_mm')==1 )+1};
        diffusion=1;
        sedimentaion=1;
        inertial=1;
        bend=1;
        contraction=1;
        tubing=1;
    else
        disp('Error: The orifice diameter in mm (d0_mm) is the inner diameter of the tube opening,the frist value of dt should equal to d0_mm')
        return
    end
else
    tubing=0;
    diffusion=0;
    sedimentaion=0;
    inertial=0;
    bend=0;
    contraction=0;
end
if(isempty(find(strcmpi(varargin,'overall_tubing'), 1))==0)
    if (isempty(find(strcmpi(varargin,'Qt'), 1))==0  &&  isempty(find(strcmpi(varargin,'L'), 1))==0)
        diffusion=1;
        Qt=varargin{find(strcmpi(varargin,'Qt')==1 )+1};
        L=varargin{find(strcmpi(varargin,'L')==1 )+1};
    else
        disp('Error:Not enough inputs for calculating the diffusion effeciency in tubing sections,need Qt or L')
        return
    end
else
    diffusion=0;
    tubing=0;
end
if(isempty(find(strcmpi(varargin,'overall_tubing'), 1))==0)
    if (isempty(find(strcmpi(varargin,'Qt'), 1))==0  &&  isempty(find(strcmpi(varargin,'L'), 1))==0 && isempty(find(strcmpi(varargin,'theta_i'), 1))==0)
        sedimentation=1;
        Qt=varargin{find(strcmpi(varargin,'Qt')==1 )+1};
        L=varargin{find(strcmpi(varargin,'L')==1 )+1};
        theta_i=varargin{find(strcmpi(varargin,'theta_i')==1 )+1};
    else
        disp('Error:Not enough inputs for calculating the sedimentation effeciency in tubing sections,need Qt or L or theta_i')
        return
    end
else
    sedimentation=0;
    tubing=0;
end
if(isempty(find(strcmpi(varargin,'overall_tubing'), 1))==0)
    if (isempty(find(strcmpi(varargin,'Qt'), 1))==0  &&  isempty(find(strcmpi(varargin,'L'), 1))==0 )
        inertial=1;
        Qt=varargin{find(strcmpi(varargin,'Qt')==1 )+1};
        L=varargin{find(strcmpi(varargin,'L')==1 )+1};
    else
        disp('Error:Not enough inputs for calculating the inertial effeciency in tubing sections,need Qt or L')
        return
    end
else
    inertial=0;
    tubing=0;
end
if isempty(find(strcmpi(varargin,'theta_Kr'), 1))==0
    bend=1;
    Qt=varargin{find(strcmpi(varargin,'Qt')==1 )+1};
    L=varargin{find(strcmpi(varargin,'L')==1 )+1};
    theta_Kr=varargin{find(strcmpi(varargin,'theta_Kr')==1 )+1};
else
    theta_Kr=0;
    bend=0;
    tubing=0;
end
if isempty(find(strcmpi(varargin,'theta_cont'), 1))==0
    contraction=1;
    Qt=varargin{find(strcmpi(varargin,'Qt')==1 )+1};
    L=varargin{find(strcmpi(varargin,'L')==1 )+1};
    theta_cont=varargin{find(strcmpi(varargin,'theta_cont')==1 )+1};
else
    theta_cont=0;
    contraction=0;
    tubing=0;
end
if tubing==1
%parameters for tubing sections
Nt=length(dt_mm)-1;%the number of the tubing section, Nt should always = length(dt)-1
dt=dt_mm*0.001;%change mm to m
%calculation
area_aver=area_aver_c(dt,Nt);%calculate average area for each tubing section aver_dt
aver_dt=aver_dt_c(area_aver,Nt);%calculate the average inner diameter of tubing
Ut=Ut_c(Qt,area_aver,Nt);%calculate average velocity in each tubing section(m/s)
%calculate the Reynolds Flow Number 
Re_tf=rho_f.*Ut.*aver_dt./K.mu;%the flow Re in the tubing 

Stkt=Stkt_c(Cc,Dp_aero,Ut,rho_p,aver_dt,K);%calculate the stokes number in each section
ef_diffusion=ef_diffusion_c(Re_tf,D,L,Qt,Sc,diffusion,Nt);%calculate the efficiency of diffusion
ef_grav=ef_grav_c(Re_tf,Vts,Ut,L,aver_dt,theta_i,sedimentation,Nt);%calculate the efficiency of sedimentation
ef_inert=ef_inert_c(Re_tf,Stkt,Ut,L,aver_dt,Qt,inertial,Nt);%calculate the efficiency of inertial
ef_inert_bend=ef_inert_bend_c(Re_tf,Stkt,theta_Kr,bend,Nt);%calculate the efficiency of inertial bend
ef_inert_cont=ef_inert_cont_c(dt,Stkt,theta_cont,contraction,Nt);%calculate the efficiency of inertial contraction
cont_at=find(theta_cont);%find which sections consider the contraction
else
ef_diffusion=1;
ef_grav=1;
ef_inert=1;
ef_inert_bend=1;
ef_inert_cont=1;
cont_at=1;
Re_tf=0;
Ut=0;
dt_mm=0;
Qt=0;
L=0;
theta_i=0;
theta_Kr=0;
theta_cont=0;
Stkt='empty';
aver_dt='empty';
end
%calculate the tubing overall efficiency
if isvector(ef_inert_cont)==1
ef_tubing=prod(ef_diffusion).*prod(ef_grav).*prod(ef_inert).*prod(ef_inert_bend).*ef_inert_cont;
else
ef_tubing=prod(ef_diffusion).*prod(ef_grav).*prod(ef_inert).*prod(ef_inert_bend).*prod(ef_inert_cont);
end
%calculate the efficiency of instruments 
%the input variables
efinstrument=fopen('instrument.txt');
instrument_scan=textscan(efinstrument,'%s %f','delimiter','=');
fclose(efinstrument);
instrument_variables={'Qa','effL_cm'};
for i_check_ins=1:length(instrument_variables)
ind_ins(i_check_ins)=isempty(find(strcmpi(varargin,instrument_variables{i_check_ins})==1,1));
end

ind_update_ins=find(ind_ins==0);

for i_up=1:length(instrument_variables)
instrument_inputs(i_up)=instrument_scan{2}(strcmpi(instrument_scan{1},instrument_variables{i_up}));
end

for i_up=1:length(ind_update_ins)
ind_var_ins(i_up)=  find(strcmpi(varargin,instrument_variables{ind_update_ins(i_up)})==1);
instrument_inputs(ind_update_ins(i_up))=varargin{ind_var_ins(i_up)+1};
end

Qa=instrument_inputs(strcmpi(instrument_variables,'Qa'));
effL_cm=instrument_inputs(strcmpi(instrument_variables,'effL_cm'));

if(4-(isempty(find(strcmpi(varargin,'DMA')==1, 1))+ isempty(find(strcmpi(varargin,'TSI_3010')==1, 1))...
        +isempty(find(strcmpi(varargin,'TSI_3025')==1, 1))+isempty(find(strcmpi(varargin,'TSI_3786')==1, 1))...
        ) >=2)
    disp('Error: can not specify instrument type, type in one instrument each time')
    return
elseif isempty(find(strcmpi(varargin,'DMA')==1, 1))==0 
    [ef_dma,mu_dma]=ef_dma_c(D,Qa,effL_cm);
    ef_instrument=ef_dma;
%     semilogx(Dp_nm,ef_dma)%for test
elseif isempty(find(strcmpi(varargin,'TSI_3010')==1, 1))==0 
    [cpc_eff]=cpc_efficiency(Dp_nm,'TSI_3010');
    ef_instrument=cpc_eff;
%     semilogx(Dp_nm,ef_instrument)%for test
elseif isempty(find(strcmpi(varargin,'TSI_3025')==1, 1))==0 
    [cpc_eff]=cpc_efficiency(Dp_nm,'TSI_3025');
    ef_instrument=cpc_eff;
%     semilogx(Dp_nm,ef_instrument)%for test
elseif isempty(find(strcmpi(varargin,'TSI_3786')==1, 1))==0 
    [cpc_eff]=cpc_efficiency(Dp_nm,'TSI_3786');
    ef_instrument=cpc_eff;
%     semilogx(Dp_nm,ef_instrument)%for test
else
    ef_instrument=1;
end


%calculate the overall inlet efficiency
ef_overall=ef_sampling.*ef_tubing.*ef_instrument;
%end of the effeciency analysis

if length(Dp_nm)>=3
semilogx(Dp_nm,ef_sampling,':b',Dp_nm,ef_tubing,'--g',Dp_nm,ef_instrument,'-.r',Dp_nm,ef_overall,'-k')
title('efficiency vs particle diameter')
ylabel('efficiency')
xlabel('Dp(nm)')
%  text(.5, .4, blue dot line is the efficiency of sampling probe)
%  text(.5, .3, green dash line is the efficiency of tubing section)
%  text(.5, .2, red dot dash line is the efficiency of instrument)
%  text(.5, .1, real black line is the overall efficiency)
end
% figure
% semilogx(Stk,ef_asp)
% loss=(1-ef_overall)*100
% loglog(Dp_nm,(1-ef_overall)*100)
% end

Input=struct('T',K.T,'P',K.P,'rho_p',rho_p,'rho_0',rho_0,'rho_in',rho_f,'charge',charge,'mu',K.mu,'mfp',K.mean_fp,...
    'U0',U0,'Q',Q,'d0_mm',d0_mm,'dout_mm',dout_mm,'dbb_mm',dbb_mm,'L_inlet_mm',L_inlet_mm...
    ,'SOT',SOT,'theta_s',theta_s,'theta_v',theta_v...
    ,'Qa',Qa,'effL_mm',effL_cm);
Properties=struct('Dp_nm',Dp_nm,'Dp_aero',Dp_aero,'Zp',Zp,'Vts',Vts,'tau_p',tau_p,'Cc',Cc,'U',U,'D',D,'Stk',Stk,'Stkt',Stkt,'Kn',Kn,'Sc',Sc,...
    'ef_aspiration',ef_asp,'ef_transmission',ef_trans,...
    'ef_diffusion',ef_diffusion,'ef_sedimentation',ef_grav,'ef_inertial',ef_inert,'ef_inertial_bend',ef_inert_bend,'ef_inertial_contraction',ef_inert_cont,'contraction_at_section',cont_at,...
    'ef_sampling_probe',ef_sampling,'ef_tubing_transport',ef_tubing,'ef_instrument',ef_instrument,'ef_overall',ef_overall);
Others=struct('Re_sampling',Re,'Re_tubing',Re_tf,'Ut',Ut,'dt_mm',dt_mm,'aver_dt',aver_dt,'Qt',Qt,'L',L,'theta_i',theta_i,'theta_Kr',theta_Kr,'theta_cont',theta_cont);
end
