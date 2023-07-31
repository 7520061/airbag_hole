function [obj,landing,pnlt] = ObjFunc_sim(para,input)

%% Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    ↓state0の初期値↓                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 数値解析に用いられる変数=微分方程式で微分値にしたいやつ
% -----------------------------------
% x(1) : エアバッグの高さ
% x(2) : 降下速度
% x(3) : 内圧
% x(4) : ガス質量
% x(5) : ガス温度
% x(6) : バッグ体積
% -----------------------------------
h0     = input.h0;      % [m]
Vc0    = input.Vc0;     % [m/s]
m0     = input.gas_m0;
T0     = input.T0;      % [m]
V0     = input.bag_V0;  % [m]

t0     = input.t0;      % Initial time [s]
state0 = [h0 Vc0 m0 T0 V0];  % initial state


switch input.sim_mode
    case 'SL'
        airbag_simulation = @AirbagSim_SL;

        [t_result, result, Gmax, Vmin] = airbag_simulation(@LandingDynamics, input, t0, state0, para);

        landing.result   = result;
        landing.t_result = t_result;
        landing.Gmax_nominal  = Gmax;
        landing.Vmin_nominal  = Vmin;

        % calculate objective
        pnlt = NaN;
        obj  = NaN;

        % Save result
        fname = [input.DIRECTORY '/data_ind'];
        n = length(dir(fname));
        save([fname '/' num2str(n-1) '_D=' num2str(para(1)) '[m]_Gmax=' num2str(Gmax) '.mat'],'landing','para','input');
        fprintf('-----------------------------------------------------------------------------------------\n')
        fprintf('D = %.2f [mm]   GMAX = %.2f [G]   VMIN = %.2f [m/s]\n' ,para(1)*1000 ,Gmax, Vmin)
        fprintf('-----------------------------------------------------------------------------------------\n')


    case 'sim'
        para = Para_sort(para);  % 距離が近い順に設計変数をソート
        airbag_simulation = @AirbagSim;

        [t_result, result, Gmax, Vmin] = airbag_simulation(@LandingDynamics, input, t0, state0, para);

        landing.result   = result;
        landing.t_result = t_result;
        landing.Gmax_nominal  = Gmax;
        landing.Vmin_nominal  = Vmin;

        % calculate objective
        pnlt = NaN;
        obj  = NaN;

        % Save result
        fname = input.DIRECTORY;
        dirname  = dir([fname '/iter*']);
        dirnum = length(dirname);
        n = length(dir([fname '/iter=' num2str(dirnum-1)]));
        save([fname '/' 'Gmax = ' num2str(Gmax) '_' num2str(n-1) '.mat'],'landing','para');
        fprintf('-----------------------------------------------------------------------------------------\n')
        fprintf('GMAX_nominal = %.2f          \t VMIN_nominal = %.2f\n' ,Gmax, Vmin)
        fprintf('-----------------------------------------------------------------------------------------\n')

end

end