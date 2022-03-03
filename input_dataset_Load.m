%% Function for making Input Data Set for learning
 
% Copyright (c) by Parikshit Pareek and Hung Nguyen, NTU Singapore 
% Date 20 Apr 2020

% Load Bus Uncertainty only
% TU=1 := Percentage uncertainty
% TU=2 := Absolute MW uncertainty


function [data_up,xx,xs,xlimit,D,rbus]=input_dataset_Load(data,Tu,N,nt,pq,sfP,sfQ,rbus,bus_DER,bus_ESS)
% lbus = data.bus(data.bus(pqbus,3)> 0,1);
% sfP=0.1; sfQ=0.1;
data_up=cell(1,N);

if Tu ==1
%     rbus=pqbus; %random bus: Bus where random power demand will be given

    pmin = data.bus(rbus,3)-abs(data.bus(rbus,3))*sfP;
    pmax = data.bus(rbus,3)+abs(data.bus(rbus,3))*sfP;

    qmin = data.bus(rbus,4)-abs(data.bus(rbus,4))*sfQ;
    qmax = data.bus(rbus,4)+abs(data.bus(rbus,4))*sfQ;
    
    % Constructing DER Bus space higher than the normal uncertainty 
    pmin(bus_DER) = pmin(bus_DER)-1; % 33-Bus total injection 3.715MW: We take 4-buses and thus will have 0.8MW capacity with 100% penetration
    % The solar is considered without providing any reactive power.
    pmin(bus_ESS) = pmin(bus_ESS)-0.8;  % ESS of 0.8MW
    pmax(bus_ESS) = pmax(bus_ESS)+0.8;  % ESS of 0.8MW
    qmin(bus_ESS) = qmin(bus_ESS)-0.9;  % Inverter of 1.2MVAR so max Q = 0.8944
    qmax(bus_ESS) = qmax(bus_ESS)-0.9;  % Inverter of 1.2MVAR
    
    % Generating training points randomly over the grid
    xlimit = [pmin' qmin';pmax' qmax'];
    % For no reverse power flow 
%     zi = xlimit <=0;
%     xlimit(zi)=0;
    D=size(xlimit,2);
    xx = rand_sample_x(N, D, xlimit); % [P's Q's]
    
  xx(end-1:end,:)=xlimit;
  xx(end-3,:)= [data.bus(rbus,3)' data.bus(rbus,4)'];
    
    % Obtaining the Q load based on the constrant power factor
    if pq==1
    pf=atan(data.bus(rbus,4)./data.bus(rbus,3));
    xx(:,D/2+1:end)=xx(:,1:D/2).*repmat(pf',[N,1]);
    end
    
    for j=1:N
        data.bus(rbus,3:4)= [xx(j,1:D/2)' xx(j,D/2+1:end)'];
        data_up{j}=data;
    end
    
%------------       Testing Sample Generation      -----------------------------
    xs = rand_sample_x(nt, D, xlimit); 
    if pq==1
        xs(:,D/2+1:end)= xs(:,1:D/2).*repmat(pf',[nt,1]);
    end
  xs(end-1:end,:)=xlimit;
  xs(end-3,:)= [data.bus(rbus,3)' data.bus(rbus,4)'];
    
        
elseif Tu==2
%      rbus=pqbus; %random bus: Bus where random power demand will be given
     
     pmin = data.bus(rbus,3)-sfP;
     pmax = data.bus(rbus,3)+sfP;
     qmin = data.bus(rbus,4)-sfQ;
     qmax = data.bus(rbus,4)+sfQ;
     

    % Generating training points randomly over the grid
     xlimit = [pmin' qmin';pmax' qmax'];
     D=size(xlimit,2);
     xx = rand_sample_x(N, D, xlimit); % [P's Q's]
    
    % Obtaining the Q load based on the constrant power factor
    if pq==1
         pf=atan(data.bus(rbus,4)./data.bus(rbus,3));
         xx(:,D/2+1:end)=xx(:,1:D/2).*repmat(pf',[N,1]);
    end
    
     
    for j=1:N
        data.bus(rbus,3:4)= [xx(j,1:D/2)' xx(j,D/2+1:end)'];
        data_up{j}=data;
    end
     
    xs = rand_sample_x(nt, D, xlimit); 
    if pq==1
        xs(:,D/2+1:end)=xs(:,1:D/2).*repmat(pf',[nt,1]);
    end

end

end