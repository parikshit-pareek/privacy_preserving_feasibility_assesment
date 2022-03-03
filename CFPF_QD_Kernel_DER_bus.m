%% Code CFPF Learning To be used for PESGM2022 Paper on PHC Analysis

% Paper: 

% Copyright: Parikshit Pareek, NTU Singapore, pare0001@ntu.edu.sg

% High Penetration at DER buses has been presented 

% data = case33bw; % The system under consideration
% 
% 
% N_training = 60; % number of samples in the tarining lattice structure 
% N_testing = 10000; % number of samples for testing

function [D_learning,Lout,Base,Mcs]  = CFPF_QD_Kernel_DER_bus(data,N_training,N_testing,sfP,sfQ,rbus,bus_DER,bus_ESS)
scale=1;
pqbus = data.bus(data.bus(:,2)==1,1); %% PQ Buses where voltage is unknown
pvbus = data.bus(data.bus(:,2)==2,1); %% PV Buses
%% --------------------------- Type of uncertainty ----------------------------------------------
%  
% All bus with nonzero load percentage Uncertainty 
% # : >  PV+PQ buses have uncertain Injection of P and  
% # : > PQ buses have uncertain Injection of Q
% Uncertainty level
% sfP=0.5; % Real power percentage fracation of base-load
% sfQ=sfP; % Reactive power percentage fracation of base-load

% kr = 2 given in code for Quadratic kernel always

[D_learning,Lout,Base]=IPF_learning_polykernel(data,pvbus,pqbus,N_training,N_testing,2,sfP,sfQ,scale,rbus,bus_DER,bus_ESS);

[~,Mcs]=MCS_output(data,D_learning.xs,D_learning.rbus,D_learning.D,N_testing,Lout,scale,pqbus);


%% Error Analysis
% Voltage Magnitude
Mcs.erV=abs(Mcs.V-Lout.muV);
Mcs.maeV = sum(abs(Mcs.erV))/N_testing;
Mcs.max_maeV = max(Mcs.maeV);
% Voltage Angle
Mcs.erTh=abs(Mcs.Thac-Lout.muTh);
Mcs.maeTh = sum(abs(Mcs.erTh))/N_testing;
Mcs.max_maeTh = max(Mcs.maeTh);



% ----------------------------------------------------------------------------------------------

function [D_learning,Lout,Base]=IPF_learning_polykernel(data,~,pqbus,N,nt,kr,sfP,sfQ,scale,rbus,bus_DER,bus_ESS)
%% Function for learning at different kernel's 
% For learning all PQ+PV buses theta and PQ bus |V| 

% Date 20 Apr 2020: Modified on 26 Apr
% Parikshit Pareek @NTUSg
% ------------------------------------------------------------------------------------------------
% ==============       Output     =============================================================
%       Lout := All output related to the learning of voltage maginitude function
% D_learning := Learning input dataset
%       Base := Base case results
%       Bopt := Optimal Bound Results
% ResControl := Voltage Control Results
pvqbus = data.bus(data.bus(:,2)<3,1);
% nbus = length(data.bus(:,1));

% pqbus = data.bus(data.bus(:,2)==1,1); %% PQ Buses where voltage is unknown

%%-     Construction of input set load uncertainty input dimension
%   -- Only load uncertainty without the powerfactor constraint  -----------
% [data_up,xx,xs,xlimit,D,rbus]=input_dataset_NetInjection(data,Tu,N,nt,2,sfP,sfQ,pqbus,pvbus);
[data_up,xx,xs,xlimit,D,rbus]=input_dataset_Load(data,1,N,nt,2,sfP,sfQ,rbus,bus_DER,bus_ESS);



%------------------ Generating Target for Training Input xx ----------------------------------------
% D_learning = Learning Dataset
parfor i=1:N
[V(:,:,i),Sg(:,:,i),Sij(:,:,i),sol{i},J{i},Y(:,:,i)]=Sampling_Jaco(data_up{i});
end

D_learning.V=V; D_learning.Sg=Sg; D_learning.Sij=Sij; 
D_learning.xx=xx; D_learning.xs=xs; D_learning.xlimit=xlimit; D_learning.D=D; % D_learning.sol=sol; 
  D_learning.data_up=data_up; D_learning.rbus=rbus;
%D_learning.Y=Y; D_learning.J=J;

%% Specify the mean, covariance and likelihood functions 
meanfunc = [];                                                  % meanZero
likfunc = @likGauss;                                            % Gaussian likelihood 
if kr == 1
    covfunc = {@covPoly,1}; % Linear (Polynomial of degree 1)  
elseif kr == 2
    covfunc = {@covPoly,2}; % Quadratic kernel (Polynomial of degree 2)  
elseif kr == 3
        covfunc = @covRQiso;  % Squared Exponental covariance function @covSEiso
elseif kr == 4
        covfunc = @covSEiso;  % Squared Exponental covariance function @covSEiso
elseif kr == 5
    covfunc = {@covMaterniso,3}; % Polynomial of degree 5
end
%% ======= Using parallel programming for learning Voltage magnitude ==========
for k=1:length(pqbus)
    yy=reshape(V(pqbus(k),1,:),[N,1])*scale;
    ytV(:,k)=yy;
        if kr < 4
            hyp = struct('mean', [], 'cov', [1; -1;1], 'lik', -1);
        else
            hyp = struct('mean',[], 'cov', [5; -1], 'lik', -1);
        end
    hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, xx, yy);
    hyp_store{k}=hyp;
    sf_lV(:,k)=exp(hyp.cov);
    [muV(:,k),s2V(:,k),~,~,~,post{k}] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, xx, yy,xs);
    alphaV(:,k)=post{k}.alpha;
end

% Defining a structure Lopt: Learning Output 
Lout.alphaV =alphaV;  Lout.muV = muV; % Lout.Smax_V=max_V;  Lout.Smin_V=min_V;  
  Lout.s2V = s2V; Lout.sf_lV=sf_lV;  Lout.ytV = ytV;  % Lout.btaV=btaV; Lout.rkhs_sqV=rkhs_sqV; 
Lout.hyps = hyp_store;  Lout.post = post;


%% ======= Using parallel programming for learning Voltage solution ==========
parfor k=1:length(pvqbus)
    yy=reshape(V(pvqbus(k),2,:),[N,1]);
    ytTh(:,k)=yy;
        if kr < 4
            hyp = struct('mean', [], 'cov', [1; -1;1], 'lik', -1);
        else
            hyp = struct('mean',[], 'cov', [5; -1], 'lik', -1);
        end
    hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, xx, yy);
%     hyp_store{k}=hyp;
    sf_lTh(:,k)=exp(hyp.cov);
    [muTh(:,k),s2Th(:,k),~,~,~,post] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, xx, yy,xs);

    alphaTh(:,k)=post.alpha;
end

% Updating a structure Lopt: Learning Output 
Lout.alphaTh =alphaTh; Lout.muTh = muTh; %Lout.btaTh=btaTh;   Lout.Smax_Th=max_Th;  Lout.Smin_Th=min_Th;  
  Lout.s2Th = s2Th; Lout.sf_lTh=sf_lTh;  Lout.ytTh = ytTh; %Lout.Kt=Kt; Lout.rkhs_sqTh=rkhs_sqTh; 



[Base.V,Base.Sg,Base.Sij,Base.sol,Base.J]=Sampling_Jaco(data);

end



%% MCS Output and error analysis for |V| at PQ buses and Theta at PV and PQ bus
function [data_upxs,Mcs]=MCS_output(data,xs,rbus,D,nt,Lout,scale,pqbus)
pvqbus = data.bus(data.bus(:,2)<3,1);
data_upxs=cell(1,nt);
scaleth=1;
% nbus = length(data.bus(:,1));

for j=1:nt
        data.bus(rbus,3:4)= [xs(j,1:D/2)' xs(j,D/2+1:end)'];
        data_upxs{j}=data;
end

 parfor i=1:nt
      [Vmc(:,:,i),Sgmc{i},Sijmc{i},solmc{i},J{i}]=Sampling_Jaco(data_upxs{i});
 end
 
 %% Error Analysis 
parfor k=1:length(pqbus)
 Vac(:,k)=Vmc(pqbus(k),1,:);
end
parfor k=1:length(pvqbus)
     Thac(:,k)=Vmc(pvqbus(k),2,:);
end
%% Error analysis of V
 Mcs.V=Vac*scale;
 Vac=Vac*scale;
 Mcs.erV_par=(abs(Vac-Lout.muV)./Vac)*100;
 
Mcs.erV_L1 = (norm(Vac*scale-Lout.muV,1)/norm(Vac*scale,1))*100;
Mcs.erV_L2 = (norm(Vac*scale-Lout.muV,2)/norm(Vac*scale,2))*100;
Mcs.erV_Linf = (norm(Vac*scale-Lout.muV,inf)/norm(Vac*scale,inf))*100;

Mcs.J = J;

%% Error in Theta 
Mcs.Thac=Thac*scaleth;
Thac=Thac*scaleth;
Mcs.erTh_par=((abs(Thac)*scaleth-abs(Lout.muTh))./(abs(Thac)))*100;
 
Mcs.erTh_L1 = (norm(abs(Thac)*scaleth-abs(Lout.muTh),1)/norm(abs(Thac)*scaleth,1))*100;
Mcs.erTh_L2 = (norm(abs(Thac)*scaleth-abs(Lout.muTh),2)/norm(abs(Thac)*scaleth,2))*100;
Mcs.erTh_Linf = (norm(abs(Thac)*scaleth-abs(Lout.muTh),inf)/norm(abs(Thac)*scaleth,inf))*100;




end









end