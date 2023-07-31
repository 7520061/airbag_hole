function  [t_result, result, Gmax, Vmin] = AirbagSim_SL(DynFun, input, t0, state0, vent_abs_para)

% -----------------------
% test
% para = rand(1,12);
% para(5) = 0;
% para = Para_sort(para);
% -----------------------

global output

state = state0;
input_temp = input;
vent_num =input.vent_num;

%% 座標から距離に変換
%中心となる穴の座標
center_abs_coord = vent_abs_para(vent_num + input.vent_close_num);
% center_abs_coord = input.close_height;

%最初に閉じる穴位置を原点とした相対座標
vent_rel_coord = (vent_abs_para(vent_num+1:end) - center_abs_coord);

% 最初に閉じる穴の中心からの距離に変換
vent_rel_D = vent_abs_para(1:vent_num);
vent_rel_H = abs(vent_rel_coord);

input_temp.s0 = vent_rel_D.^2 / 4 *pi; %各穴面積
input_temp.n0 = input.n0_ind * ones(vent_num,1); %穴数





%% dynamics solver
% 次のresultがグラフ書くため用の行列で，1行目にとりあえず初期値を保存．次からはi行目にiのときの結果を保存される(f_landingのときに実行)．
% グラフ化したい結果を設定
% -----------------------------------
% result(1) : エアバッグの高さ
% result(2) : 降下速度
% result(3) : 内圧(だるいけkPaで表示)
% result(4) : ガス質量
% result(5) : ガス温度
% result(6) : バッグ体積
% result(7) : 荷重倍数G
% result(8) : 接触面積
% result(9) : マッハ数
% result(10) : ベントホール面積
% result(11) : 排気速度
% -----------------------------------
result(1,:) = [state(1) state(2) input.p0/1000 state(3) state(4) state(5) (input.acc+input.g)/input.g input.bag_S input.mae input.vent_s input.ve input_temp.Vflag -input.d_mass];
t_result(1) = t0;      % 時間

i = 1;
Gmax = 0;
Vmin = 100;

for t = t0 : input.dt : input.t_end
    i = i+1;
    %       -----(  ・ω・)ノこっからが本命じゃー！！------
    [t_n, state_n] = ode45(@(t, state)DynFun(state, input_temp), [t t+input.dt], state);
    %       -------------------------------------------------
    % 解説
    % odeは微分方程式を解く際に使う関数
    % [T,X] = ode45(@function,tspan, x0)の形で表される
    % tspanは積分範囲のことで，何秒から何秒の区間で積分しますか？という意味
    % x0は微分方程式functionに与える初期値（今回の場合はt秒のときのxの値という意味
    % functionで微分法定式を解く⇒ TのときのXの値がわかる
    % Tの行列には積分区間をさらに細かく刻んだ時間が保存している（今回は40分割で40行くらいに保存されてる）
    % Xの行列にはそのTに対応した結果が保存されている
    % まあ最終的に積分区間の最後であるt+dtのときのxの値がわかればいいので行列の最後の値だけわかればええ
    % X(end,:)←これがほしい

    state = state_n(end,:); % ここでｘの値に更新してる　次の計算ではこの値が微分方程式に代入される
    t_result(i) = t_n(end,1);
    result(i,:) = [state(1) state(2) output.bag_p/1000 state(3) state(4) state(5) output.G output.bag_S output.Mae output.se output.ve input_temp.Vflag output.d_mass];

    input_temp.bag_H = output.bag_h;
    input_temp.bag_p = output.bag_p;
    input_temp.acc = output.acc;

    % バッグ内圧が更新されたけベントホールが開口するかどうか判断　ちゃんと一回空いたらもう閉じないような設定にしてる
    % input.vent_close_numで指定された穴から閉じる

    if (result(i,3) > input.break_pressure || result(i,10) > 0) && result(i,1) <= input.bag_Hs %円筒部に入る

        hz = (input.bag_Hs - result(i,1)) / 2; % 円筒部内で折り畳まれた距離


        if hz <= input.bag_Hs/2 - abs(center_abs_coord)  % phase1...折り畳まれている部分が穴の中心から上下に進む，上下どちらかの端に到達するまで

            cnt = input.vent_close_num;

            if hz <= vent_rel_D(input.vent_close_num)

                if input_temp.n0(cnt) == 0
                    continue
                end

                input_temp.n0(cnt) = input.n0_ind * vent_block_ratio(hz,vent_rel_H(cnt),vent_rel_D(cnt),'2');
            end

            for j = 1:1:vent_num

                if j == cnt || input_temp.n0(j) == 0
                    continue
                end

                if vent_rel_H(j) - vent_rel_D(j)/2 <= hz && hz <= vent_rel_H(j) + vent_rel_D(j)/2
                    input_temp.n0(j) = input.n0_ind * vent_block_ratio(hz,vent_rel_H(j),vent_rel_D(j),'1');
                end

                if vent_rel_H(j) + vent_rel_D(j)/2 < hz
                    input_temp.n0(j) = 0;
                end

            end


        else  % phase2...折り畳まれている部分が下方向に進む

            if input.vent_close_num > vent_num/2
                phase2_end_rel_coord = -(input.bag_Hs/2 + center_abs_coord);
                phase2_hz = phase2_end_rel_coord - hz;

            else
                phase2_end_rel_coord = input.bag_Hs/2 + center_abs_coord;
                phase2_hz = phase2_end_rel_coord + hz;

            end

            for k = 1:1:vent_num

                if input_temp.n0(k) == 0
                    continue
                end

                if vent_rel_coord(k) - vent_rel_D(k)/2 <= phase2_hz && phase2_hz <= vent_rel_coord(k) + vent_rel_D(k)/2
                    input_temp.n0(k) = input.n0_ind * vent_block_ratio(phase2_hz,vent_rel_coord(k),vent_rel_D(k),'1');
                end

                if vent_rel_coord(k) - vent_rel_D(k)/2 > phase2_hz
                    input_temp.n0(k) = 0;
                end
            end
        end

    end


    input_temp.Vflag = sum(input_temp.n0);
    input_temp.vent_s = input_temp.s0 * input_temp.n0;


    if result(i,1) > input_temp.bag_dh && Gmax <= result(i,7)
        Gmax = max(result(i,7));
    end

    if result(i,1) > input_temp.bag_dh && Vmin >= abs(result(i,2))
        Vmin = min(abs(result(i,2)));
    end

    % エアバッグがぺったんこになったら強制終了（(；゜Д゜)
    if result(i,1) <= input_temp.bag_dh
        break;
    end

% 高さがマイナスになるとGが発散する為それを除いたうえでGmax,Vminを見る
% terget_res = result(result(:,1)>0,:);



end

