% %% ----------------------- 1. Learn the V = f(s) using Quadratic Kernel ------------------------------------------

function [Res,D_learning,Lout,Mcs] = battery_feasibility_space(data,N_training,N_testing,sfP,bus_DER,bus_ESS,hr_idx,sf_lb_P)
% sfP=0.3; % Real power percentage fracation of base-load
sfQ=sfP; % Reactive power percentage fracation of base-load 
% bus_DER = [18,22,29]; % DER NODES (-1 to remove slack bus when used for index)
pqbus = data.bus(data.bus(:,2)==1,1); %% PQ Buses where voltage is unknown
rbus = pqbus;  
% bus_ESS = [15]-1; % ESS NODES (-1 to remove slack bus when used for index)

[D_learning,Lout,~,Mcs]  = CFPF_QD_Kernel_DER_bus(data,N_training,N_testing,sfP,sfQ,rbus,bus_DER,bus_ESS);

D = D_learning.D;
%% ====================== Obtaining the Voltage value at different points ===================================================
% Idea:
% -- Generate samples of Load, DER injection and control limits of reactive power
% -- Check if there exist one control sample for which the capacity is admissible
 
% --------------------- Sl : Set of Uncertain Load points between slb and sub ---------------------------------------------------------------------
% sf_lb_P = 0.1; 
sf_lb_Q = sf_lb_P; 
slb = [ data.bus(D_learning.rbus,3)-abs(data.bus(D_learning.rbus,3))*sf_lb_P;data.bus(D_learning.rbus,4)-abs(data.bus(D_learning.rbus,4))*sf_lb_Q];
sub = [ data.bus(D_learning.rbus,3)+abs(data.bus(D_learning.rbus,3))*sf_lb_P;data.bus(D_learning.rbus,4)+abs(data.bus(D_learning.rbus,4))*sf_lb_Q];
s_limits = [slb sub];
% so = [data.bus(D_learning.rbus,3);data.bus(D_learning.rbus,4)];
eps = 0.01;
Ns = round(log(1/10^(-6))/log(1/(1-eps)));% Number of samples needed for load space
% Nq = 1000;
Ne = 100; % Number of points in each dimension delta =  (x_max-x_min)/100
for i = 1:D
    X_all(:,i) = linspace(slb(i),sub(i),Ne);
end
rxt = randi(Ne,Ns,D);
% rxs = randi(Ne,N_xs,D);
for i = 1: Ns
    for j = 1:D
        Sl(i,j) = X_all(rxt(i,j),j); % Each column of this is a random vector for training
    end
end
load('pv.mat');
% ---------------------------------------------Generating DER Sample Sets for a given level of penetration _Px_ ------------------------------------------
% ---------------------------------------- Making the control samples of reactive power at the ndoes _Qcap_ and _Qc_--------------------------------
for i = 1:1
%      total_DER = 1; %MW  sum(so(1:D/2))*(Cap_DER(i)/100); % Total DER in system as per percentage of base load
total_DER =  pv(hr_idx);
%     if hr_idx <= 7
%     total_DER = 0;
%     elseif hr_idx >= 18
%        total_DER = 0;
%     else
%       total_DER = abs(6-hr_idx);
%     end
        ind_cap = total_DER/length(bus_DER);    % Individual capacity 

%     inv_cap = B_cap*1.25; % Inverter capacity is 1.25 times of MW capacity => Q capacity is 0.75 of MW capacity
    
    
    for j = 1:length(bus_DER)
        Xp(:,j) = linspace(0,ind_cap,Ne);
%         Xpb(:,j) = linspace(-B_cap,B_cap,Ne*2);
    end
    
    rp = randi(Ne,Ns,length(bus_DER));
    
    for ii = 1: Ns
        for j = 1:length(bus_DER)
            Px(ii,j,i) = Xp(rp(ii,j),j); % Each layer of this means one DER_Capacity, each row is one sample and total 4 columns (Equal to #DER)
        end
    end
          
end

% ----------------------------- Generate points of feasibility testing for the ESS -------------------------------------------------------
D_grid = 30;
Pb = linspace(-0.8,0.8,D_grid);
Qb = linspace(-0.9,0.9,D_grid);

[XPb, XQb] = meshgrid(Pb,Qb);

%  ------------- Quadratic Form V = x'M x +Nx+c^2 --------------------------------------- 
% ==================================================================================================
alphaV = Lout.alphaV;
l=Lout.sf_lV(1,:);
xx = D_learning.xx;
% xs = D_learning.xs;
D = D_learning.D;
tau = Lout.sf_lV(3,:);
c = Lout.sf_lV(2,:);

parfor k=1:length(pqbus)
alpha_l4(:,k)= alphaV(:,k)/(l(k)^4);
end
parfor i=1:N_training
    M(:,:,i) = xx(i,:)'*xx(i,:);
end
M_alpha = zeros(D,D,length(pqbus));
for k=1:length(pqbus)
    for i=1:N_training
        M_alpha(:,:,k)=M_alpha(:,:,k) + alpha_l4(i,k)*M(:,:,i);
    end
    N_alpha(k,:) = (2*alphaV(:,k)'*(c(k))*xx)/(l(k)^2);
end

z = zeros(Ns,D/2);
zb = zeros(Ns,D/2);

% Construction of Combined Sample Batches
d = 1; % Only one value in DER capacity array
zb(:,bus_DER) = Px(:,:,d);
SP(:,:,d) = Sl-[zb z]; % 3-D matrix: [Ns x D x Ns  ]

% Calculating the voltage values for the 'SP' input where 'S' is load and 'P' is the renewable injectio 

     for k=1:length(pqbus)
        for j=1:Ns
             V_SP(j,k)= tau(k)^2*(SP(j,:,d)*M_alpha(:,:,k)*SP(j,:,d)'+N_alpha(k,:)*SP(j,:,d)'+alphaV(:,k)'*((c(k)^2)*ones(N_training,1)));
        end
    end



zpb = zeros(1,D/2);
zqb = zeros(1,D/2);
% Calculating the voltage values for the Pb-Qb: Control feasibility testing from battery point of view
for i = 1:D_grid
    for j = 1:D_grid
        zqb(:,bus_ESS) = XQb(i,j);
        zpb(:,bus_ESS) = XPb(i,j);
    for k=1:length(pqbus)
        V_B(i,j,k)= tau(k)^2*([zpb zqb]*M_alpha(:,:,k)* [zpb zqb]'+N_alpha(k,:)* [zpb zqb]');
    end
    end
end

% Calculating the final voltage value Vsp+Vq for each sample
for j = 1:D_grid
    for i  = 1:D_grid
         V_3D{i,j} =  V_SP + repmat(reshape(V_B(i,j,:),[1,32]),[Ns,1]);
    end
end
% Now each of the slice (3rd dimension) of V_3D tells us if that perticular Qc has able to provide hosting 
% We need only one Qc (one layer of V_3D) to be completely within limits of voltage 
for i = 1:D_grid
    for j = 1:D_grid
        V_min(i,j) = min(min(V_3D{i,j}'));
        V_max(i,j) = max(max(V_3D{i,j}'));
    end
end

f_min = V_min > 0.90; % Feasible minimum index 
f_max = V_max < 1.10; % Feasible maximum index 

f = (f_min + f_max) == 2;

XPb_f = XPb(f); % Feasible Real Power Injection By battery
XQb_f = XQb(f); % Feasible Reactive Power Injection By battery

% toc
% 
% hr_idx = 1;
Res.f = f;
Res.Xb_f = [XPb_f XQb_f];
% %  Calculate the oltage Control Performance Index and Price
% Res{hr_idx}.VCP = V_max(f)-V_min(f);
Res.Vmax_f = V_max(f);
Res.Vmax_f = V_min(f);
% %% Saving results
Res.s_limits = s_limits;
Res.ind_cap = ind_cap;

Res.s_limits(bus_DER,1) = s_limits(bus_DER,1)-ind_cap;
% Res{hr_idx}.Mcs = Mcs;
Res.XPb = XPb;
Res.XQb = XQb;
end
% % end
% plot(Res{hr_idx}.Xb_f(:,1),Res{hr_idx}.Xb_f(:,2),'*')
% hold on
% subplot(1,3,1)
% k = boundary(Res{hr_idx}.Xb_f(:,1),Res{hr_idx}.Xb_f(:,2));
% plot(Res{hr_idx}.Xb_f(k,1),Res{hr_idx}.Xb_f(k,2))
% hold on


