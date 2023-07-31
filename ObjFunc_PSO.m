function [obj,landing,pnlt] = ObjFunc_PSO(para,input)

%% Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    ↓state0の初期値↓                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 数値解析に用いられる変数=微分方程式で微分値にしたいやつ
% -----------------------------------
% x(1) : エアバッグの高さ
% x(2) : 降下速度
% x(3) : ガス質量
% x(4) : ガス温度
% x(5) : バッグ体積
% -----------------------------------
h0     = input.h0;      % [m]
Vc0    = input.Vc0;     % [m/s]
m0     = input.gas_m0;
T0     = input.T0;      % [m]
V0     = input.bag_V0;  % [m]

t0     = input.t0;      % Initial time [s]
state0 = [h0 Vc0 m0 T0 V0];  % initial state

para = Para_sort(para,input.vent_num);  % 距離が近い順に設計変数をソート
vent_D_sorted = para(1:input.vent_num);
vent_H_sorted = para(input.vent_num+1:end);


[t_result, result, Gmax, Vmin] = AirbagSim(@LandingDynamics, input, t0, state0, para);


landing.result   = result;
landing.t_result = t_result;
landing.Gmax_nominal  = Gmax;
landing.Vmin_nominal  = Vmin;
Gmax_threesigma = 100;
Gmax_worst = 100;

%% calculate objective
% ベントホールの有無を判定
vent_flag = [0 0 0 0 0 0];
for j = 1:1:6
    if vent_D_sorted(j) ~= 0
        vent_flag(j) = j;
    end
end

B = vent_flag == 0;
vent_H_sorted(B) = [];
vent_D_sorted(B) = [];
n = numel(vent_D_sorted);

% %ベントホール間距離
% dif = ones(n-2,1);
% pnlt_temp = ones(n-2,1);
% for j = 1:1:n-2
%     dif(j) = ((vent_H_sorted(j) - vent_D_sorted(j)) - (vent_H_sorted(j+2) + vent_D_sorted(j+2)))/2;
%     if dif(j) <= 0.010 %0.025 % 0.010
%         pnlt_temp(j) = (0.025 - dif(j)) * 1e+20;
%     elseif dif(j) < 0.025
%         pnlt_temp(j) = (0.025 - dif(j)) * 100;
%     else
%         pnlt_temp(j) = 0;
%     end
% end

% pnlt1 = sum(pnlt_temp);
% JEA1  = pnlt1;

% pnlt1 = 0;
% JEA1 = 0;

% 円筒部端からベントホールまでの距離
d_edge = zeros(n,1);
pnlt_temp = zeros(n,1);

for j =1:1:n
    d_edge(j) = abs(input.bag_Hs - (vent_H_sorted(j) + vent_D_sorted(j)));
    if d_edge(j) <= 0.060
        pnlt_temp(j) = d_edge(j) * 1e+20;
    else
        pnlt_temp(j) = 0;
    end
end

pnlt1 = sum(pnlt_temp);
JEA1  = pnlt1;

% 降下速度
if Vmin >= 1.0
    pnlt2 = 1e+15;
else
    pnlt2 = 0;
end
    
JEA2  = abs(result(end,2)) * pnlt2;

% 荷重倍数
if Gmax > input.lf_limit || pnlt1 >= 1e+10
    pnlt3 = Gmax * 1e+08;
    JEA3 =  Gmax * pnlt3;
    % ノミナルがそもそもダメなので3σまで見ない
    pnlt4 = 0;
    JEA4  = 0;
else 
    pnlt3 = 0;
    JEA3 =  0;
    %3σ荷重倍数
    [t_worst, result_worst, Gmax_worst, Gmax_threesigma, maxload, maxparam] = error_propagation(@LandingDynamics, input, t0, state0, para, Gmax);

    landing.result_worst = result_worst;
    landing.t_worst = t_worst;
    landing.Gmax_worst  = Gmax_worst;
    landing.Gmax_threesigma = Gmax_threesigma;
    landing.maxload_error   = maxload;
    landing.maxparam_error  = maxparam;
    
    if Gmax_threesigma <= input.lf_limit
        pnlt4 = 1; 
    else
        pnlt4 = 1e+04;
    end
    JEA4 =  (Gmax_threesigma - input.lf_limit) * pnlt4;
end
pnlt = [pnlt1 pnlt2 pnlt3 pnlt4];
obj  = JEA1 + JEA2 + JEA3 + JEA4;

%% Save result
fname = input.DIRECTORY;
dirname  = dir([fname '/iter*']);
dirnum = length(dirname);
fprintf('iter = %.f  ID = %.f   GMAX_nominal = %.2f          \t VMIN_nominal = %.2f   obj = %.f\n' ,dirnum-1 ,get(getCurrentTask(),'ID'), Gmax, Vmin, obj)

if Gmax <= input.lf_limit && Vmin <= 1.0% dirnum<=input.nGen && Gmax_max <= input.lf_limit
    n = length(dir([fname '/iter=' num2str(dirnum-1)]));
    save([fname '/iter=' num2str(dirnum-1) '/obj = ' num2str(obj) '_' 'Gmax = ' num2str(Gmax) '_' '3sigma = ' num2str(Gmax_threesigma) '_' num2str(get(getCurrentTask(),'ID')) '_'  num2str(n-1) '.mat'],'obj','landing','pnlt','para');
    fprintf('-----------------------------------------------------------------------------------------\n')
    fprintf('iter = %.f  ID = %.f   GMAX_nominal = %.2f          \t VMIN_nominal = %.2f   obj = %.f\n' ,dirnum-1 ,get(getCurrentTask(),'ID'), Gmax, Vmin, obj)
    fprintf('GMAX_worst = %.2f  GMAX_3sigma = %.2f\n',Gmax_worst,Gmax_threesigma)
    fprintf('-----------------------------------------------------------------------------------------\n')
end

end