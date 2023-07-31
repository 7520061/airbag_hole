function [obj,landing,pnlt] = ObjFunc_PSO2(para,input)

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


switch input.calc_mode
    case 1
        %% sort
        para = Para_sort(para);  % 距離が近い順に設計変数をソート
        para_num = length(para);
        vent_D_sorted = para(1:para_num/2);
        vent_H_sorted = para(para_num/2+1:end);

        %% airbag simulation
        % 様々なつぶれ方で評価

        % method 1
        input.vent_close_method = input.vent_close_method1;
        [t_result_1, result_1, Gmax_1, Vmin_1] = AirbagSim(@LandingDynamics, input, t0, state0, para);

        landing.method1.result   = result_1;
        landing.method1.t_result = t_result_1;
        landing.method1.Gmax_nominal  = Gmax_1;
        landing.method1.Vmin_nominal  = Vmin_1;

        % method 2
        input.vent_close_method = input.vent_close_method2;
        [t_result_2, result_2, Gmax_2, Vmin_2] = AirbagSim(@LandingDynamics, input, t0, state0, para);

        landing.method2.result   = result_2;
        landing.method2.t_result = t_result_2;
        landing.method2.Gmax_nominal  = Gmax_2;
        landing.method2.Vmin_nominal  = Vmin_2;

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

        % JEA1 円筒部端からベントホールまでの距離
        d_edge = zeros(n,1);
        pnlt1_temp = zeros(n,1);

        for j =1:1:n
            d_edge(j) = abs(input.bag_Hs - (vent_H_sorted(j) + vent_D_sorted(j)));
            if d_edge(j) <= 0.060
                pnlt1_temp(j) = d_edge(j) * 1e+20;
            else
                pnlt1_temp(j) = 0;
            end
        end

        pnlt1 = sum(pnlt1_temp);
        JEA1  = pnlt1;


        % JEA2 降下速度
        % 降下速度を違反している数だけペナルティを与える
        Vmin = [Vmin_1,Vmin_2];

        pnlt2_temp = Vmin >= 1.0; %Vminが1.0以上でpnlt2_temp=1,未満でpnlt2_temp=0

        pnlt2 = pnlt2_temp .* 1e+15; %違反してたら1e+15,してなかったら0
        JEA2  = sum(abs(Vmin) .* pnlt2);


        % JEA3 荷重倍数
        Gmax = [Gmax_1,Gmax_2];

        pnlt3_temp = Gmax > input.lf_limit; %Gmaxが7より上でpnlt3_temp=1,以下でpnlt3_temp=0

        pnlt3 = pnlt3_temp .* 1e+08;
        JEA3 =  sum(Gmax .* pnlt3);


        % JEA4 3σ荷重倍数(未実装)
        %
        % [t_worst, result_worst, Gmax_worst, Gmax_threesigma, maxload, maxparam] = error_propagation(@LandingDynamics, input, t0, state0, para, Gmax);
        %
        % landing.result_worst = result_worst;
        % landing.t_worst = t_worst;
        % landing.Gmax_worst  = Gmax_worst;
        % landing.Gmax_threesigma = Gmax_threesigma;
        % landing.maxload_error   = maxload;
        % landing.maxparam_error  = maxparam;
        %
        % if Gmax_threesigma <= input.lf_limit
        %     pnlt4 = 1;
        % else
        %     pnlt4 = 1e+04;
        % end
        % JEA4 =  (Gmax_threesigma - input.lf_limit) * pnlt4;
        pnlt4 = 0;
        JEA4 = 0;

        pnlt = [pnlt1 pnlt2 pnlt3 pnlt4];
        obj  = JEA1 + JEA2 + JEA3 + JEA4;

        %% Save result
        fname = input.DIRECTORY;
        dirname  = dir([fname '/iter*']);
        dirnum = length(dirname);
        fprintf('iter = %.f  ID = %.f   Gmax_1 = %.2f   Gmax_2 = %.2f          \t Vmin_1 = %.2f  Vmin_2 = %.2f   obj = %.f\n' ,dirnum-1 ,get(getCurrentTask(),'ID'), Gmax_1, Gmax_2, Vmin_1, Vmin_2, obj)


        if  sum(pnlt2) == 0 && sum(pnlt3) == 0 %7G以下かつ降下速度1m/s以下に抑えられている個体を表示
            n = length(dir([fname '/iter=' num2str(dirnum-1)]));
            save([fname '/iter=' num2str(dirnum-1) '/obj = ' num2str(obj) '_' 'Gmax_1 = ' num2str(Gmax_1) '_' 'Gmax_2 = ' num2str(Gmax_2) '_' num2str(get(getCurrentTask(),'ID')) '_'  num2str(n-1) '.mat'],'obj','landing','pnlt','para');
            fprintf('---------------------------------------------------------------------------------------------------------\n')
            fprintf('iter = %.f  ID = %.f   Gmax_1 = %.2f   Gmax_2 = %.2f          \t Vmin_1 = %.2f  Vmin_2 = %.2f   obj = %.f\n' ,dirnum-1 ,get(getCurrentTask(),'ID'), Gmax_1, Gmax_2, Vmin_1, Vmin_2, obj)
            fprintf('---------------------------------------------------------------------------------------------------------\n')
        end
    case 2
        %% 上半分で対称になるように配置
        des_num = length(para); %設計変数の数

        vent_H = para(des_num/2+1:end);
        vent_H = [vent_H,flip(input.bag_Hs-vent_H)];

        vent_D = para(1:des_num/2);
        vent_D = [vent_D,flip(vent_D)];

        para_calc = [vent_D,vent_H];

        %% sort
        para_calc = Para_sort(para_calc);  % 距離が近い順に設計変数をソート
        para_calc_num = length(para_calc);
        vent_D_sorted = para_calc(1:para_calc_num/2);
        vent_H_sorted = para_calc(para_calc_num/2+1:end);


        %% airbag simulation
        input.vent_close_method = '中腹から順次閉じる_面積比から閉じ量求める';%中腹でも上下でも同じ
        [t_result, result, Gmax, Vmin] = AirbagSim(@LandingDynamics, input, t0, state0, para_calc);

        landing.result   = result;
        landing.t_result = t_result;
        landing.Gmax_nominal  = Gmax;
        landing.Vmin_nominal  = Vmin;
        Gmax_threesigma = 100;
        Gmax_worst = 100;

        %% calculate objective
        % ベントホールの有無を判定
        vent_flag = [0 0 0 0 0 0];
        for j = 1:1:para_calc_num/2
            if vent_D_sorted(j) ~= 0
                vent_flag(j) = j;
            end
        end

        B = vent_flag == 0;
        vent_H_sorted(B) = [];
        vent_D_sorted(B) = [];
        n = numel(vent_D_sorted);

        % 円筒部端からベントホールまでの距離
        d_edge = zeros(n,1);
        pnlt_temp = zeros(n,1);

        for j =1:1:n
            d_edge(j) = input.bag_Hs - (vent_H_sorted(j) + vent_D_sorted(j));
            if d_edge(j) <= 0.060
                pnlt_temp(j) = abs(d_edge(j)) * 1e+20;
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
            [t_worst, result_worst, Gmax_worst, Gmax_threesigma, maxload, maxparam] = error_propagation(@LandingDynamics, input, t0, state0, para_calc, Gmax);

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

end