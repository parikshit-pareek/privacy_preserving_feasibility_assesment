%% Code for recoding the video and releasing on Github
clear
clc



load('isu_data.mat')
Q_j1 = ISU_Data.Q32_june(:,:,1)'/20;
P_j1 = ISU_Data.P32_june(:,:,1)'/20;

for  hr_idx = 1:24
data = case33bw; % The system under consideration


data.bus(2:end,3) = P_j1(:,hr_idx);
data.bus(2:end,4) = Q_j1(:,hr_idx);

N_training = 400; % number of samples in the tarining lattice structure 
N_testing = 2000; % number of samples for testing

sfP=0.3; % For Training --> Fracation of base-load for GENERAL UNCERTAINTY @ All Nodes
bus_DER = [18,22,31]-1;
bus_ESS = [15]-1;

sf_lb_P = 0.1; % < sfP ALWAYS; For  Feasibility Space Idenfication --> Fracation of base-load for GENERAL UNCERTAINTY @ All Nodes
[R{hr_idx},D_learning{hr_idx},Lout{hr_idx},Mcs{hr_idx}] = battery_feasibility_space(data,N_training,N_testing,sfP,bus_DER,bus_ESS,hr_idx,sf_lb_P);

k = boundary(R{1,hr_idx}.Xb_f(:,1), R{1,hr_idx}.Xb_f(:,2));
figure
xlim([-0.85 0.85])
ylim([-0.905 0.905])
plot( R{1,hr_idx}.Xb_f(k,1), R{1,hr_idx}.Xb_f(k,2),'--')
ylabel('Reactive Power')
xlabel('Real Power')
title('Negative is Battery Discharging')
end

