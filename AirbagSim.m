function  [t_result, result, Gmax, Vmin] = AirbagSim(DynFun, input, t0, state0, para)

global output

state = state0;
input_temp = input;




for j = 1 : 1 : 6
    if para(j) == 0.000
        if j == 1
            input_temp.an0 =0;
        elseif j == 2
            input_temp.bn0 =0;
        elseif j == 3
            input_temp.cn0 =0;
        elseif j == 4
            input_temp.dn0 =0;
        elseif j == 5
            input_temp.en0 =0;
        elseif j == 6
            input_temp.fn0 =0;
        end
    end
end

input_temp.as0 = para(1)^2 / 4 * pi ;%* input_temp.an0;
input_temp.bs0 = para(2)^2 / 4 * pi ;%* input_temp.bn0;
input_temp.cs0 = para(3)^2 / 4 * pi ;%* input_temp.cn0;
input_temp.ds0 = para(4)^2 / 4 * pi ;%* input_temp.dn0;
input_temp.es0 = para(5)^2 / 4 * pi ;%* input_temp.en0;
input_temp.fs0 = para(6)^2 / 4 * pi ;%* input_temp.fn0;
% input_temp.gs0 = para(13)^2 / 4 * pi; %常に開いている穴
input_temp.gs0 = 0;

% vent_h_check = [para(7) para(8) para(9) para(10) para(11) para(12)];
% vent_h_min = min(vent_h_check);
% input.bag_dh = para(12) - para(6);    % 最小バッグ高さ[m]

input_temp.Vflag = 0;

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

    % ベントホール閉じ方
    switch input.vent_close_method
        case '順次閉じる'
            if result(i,3) > input.break_pressure || result(i,10) > 0

                if para(7) >= result(i,1)
                    input_temp.an0 = 0;   
                    input_temp.gn0 = 0;   
                end
                if para(8) >= result(i,1)
                    input_temp.bn0 = 0;   
                    input_temp.hn0 = 0;   
                end
                if para(9) >= result(i,1)
                    input_temp.cn0 = 0;   
                    input_temp.in0 = 0;   
                end
                if para(10) >= result(i,1)
                    input_temp.dn0 = 0;   
                    input_temp.jn0 = 0;   
                end
                if para(11) >= result(i,1)
                    input_temp.en0 = 0;   
                    input_temp.kn0 = 0;  
                end
                if para(12) >= result(i,1)
                    input_temp.fn0 = 0;   
                    input_temp.ln0 = 0;   
                end

                input_temp.Vflag  = input_temp.an0 + input_temp.bn0 + input_temp.cn0 + input_temp.dn0 + input_temp.en0 + input_temp.fn0; % + input_temp.gn0 + input_temp.hn0 + input_temp.in0 + input_temp.jn0 + input_temp.kn0 + input_temp.ln0;
                input_temp.vent_s = input_temp.as0*input_temp.an0 + input_temp.bs0*input_temp.bn0 + input_temp.cs0*input_temp.cn0...
                    + input_temp.ds0*input_temp.dn0 + input_temp.es0*input_temp.en0 + input_temp.fs0*input_temp.fn0;
            end
        case '中腹から順次閉じる'
            if result(i,3) > input.break_pressure || result(i,10) > 0

                if input.bag_Hs - para(12) >= result(i,1)
                    input_temp.fn0 = 0;
                end

                if input.bag_Hs - para(11) >= result(i,1)
                    input_temp.en0 = 0;
                end

                if input.bag_Hs - para(10) >= result(i,1)
                    input_temp.dn0 = 0;
                end

                if input.bag_Hs - para(9) >= result(i,1)
                    input_temp.cn0 = 0;
                end

                if input.bag_Hs - para(8) >= result(i,1)
                    input_temp.bn0 = 0;
                end

                if input.bag_Hs - para(7) >= result(i,1)
                    input_temp.an0 = 0;
                end

                
                input_temp.Vflag  = input_temp.an0 + input_temp.bn0 + input_temp.cn0 + input_temp.dn0 + input_temp.en0 + input_temp.fn0; % + input_temp.gn0 + input_temp.hn0 + input_temp.in0 + input_temp.jn0 + input_temp.kn0 + input_temp.ln0;
                input_temp.vent_s = input_temp.as0*input_temp.an0 + input_temp.bs0*input_temp.bn0 + input_temp.cs0*input_temp.cn0...
                    + input_temp.ds0*input_temp.dn0 + input_temp.es0*input_temp.en0 + input_temp.fs0*input_temp.fn0+input_temp.gs0;
            end
            
        case '中腹から順次閉じる_面積比から閉じ量求める'
            hz = input.bag_Hs - result(i,1); %エアバッグが閉じている高さ

            if (result(i,3) > input.break_pressure || result(i,10) > 0) && result(i,1) <= input.bag_Hs 

                if input.bag_Hs - para(12) + para(6) >= result(i,1) && result(i,1) > input.bag_Hs - para(12) - para(6)
                    input_temp.fn0 = input.fn0 * vent_block_ratio3(hz,para(12),para(6),'1');

                elseif input.bag_Hs - para(12) - para(6) >= result(i,1)
                    input_temp.fn0 = 0;
                end


                if input.bag_Hs - para(11) + para(5) >= result(i,1) && result(i,1) > input.bag_Hs - para(11) - para(5)
                    input_temp.en0 = input.en0 * vent_block_ratio3(hz,para(11),para(5),'1');

                elseif input.bag_Hs - para(11) - para(5) >= result(i,1)
                    input_temp.en0 = 0;
                end


                if input.bag_Hs - para(10) + para(4) >= result(i,1) && result(i,1) > input.bag_Hs - para(10) - para(4)
                    input_temp.dn0 = input.dn0 * vent_block_ratio3(hz,para(10),para(4),'1');

                elseif input.bag_Hs - para(10) - para(4) >= result(i,1)
                    input_temp.dn0 = 0;
                end


                if input.bag_Hs - para(9) + para(3) >= result(i,1) && result(i,1) > input.bag_Hs - para(9) - para(3)
                    input_temp.cn0 = input.cn0 * vent_block_ratio3(hz,para(9),para(3),'1');

                elseif input.bag_Hs - para(9) - para(3) >= result(i,1)
                    input_temp.cn0 = 0;
                end


                if input.bag_Hs - para(8) + para(2) >= result(i,1) && result(i,1) > input.bag_Hs - para(8) - para(2)
                    input_temp.bn0 = input.bn0 * vent_block_ratio3(hz,para(8),para(2),'1');

                elseif input.bag_Hs - para(8) - para(2) >= result(i,1)
                    input_temp.bn0 = 0;
                end


                if input.bag_Hs - para(7) + para(1) >= result(i,1) && result(i,1) > input.bag_Hs - para(7) - para(1)
                    input_temp.an0 = input.an0 * vent_block_ratio3(hz,para(7),para(1),'1');

                elseif input.bag_Hs - para(7) - para(1) >= result(i,1)
                    input_temp.an0 = 0;
                end

%                 input_temp.an0 = input.an0 * vent_block_ratio(input.bag_Hs,para(1),para(7) ,input_temp.as0,result(i,1));
%                 input_temp.bn0 = input.bn0 * vent_block_ratio(input.bag_Hs,para(2),para(8) ,input_temp.bs0,result(i,1));
%                 input_temp.cn0 = input.cn0 * vent_block_ratio(input.bag_Hs,para(3),para(9) ,input_temp.cs0,result(i,1));
%                 input_temp.dn0 = input.dn0 * vent_block_ratio(input.bag_Hs,para(4),para(10),input_temp.ds0,result(i,1));
%                 input_temp.en0 = input.en0 * vent_block_ratio(input.bag_Hs,para(5),para(11),input_temp.es0,result(i,1));
%                 input_temp.fn0 = input.fn0 * vent_block_ratio(input.bag_Hs,para(6),para(12),input_temp.fs0,result(i,1));

                input_temp.Vflag  = input_temp.an0 + input_temp.bn0 + input_temp.cn0 + input_temp.dn0 + input_temp.en0 + input_temp.fn0; % + input_temp.gn0 + input_temp.hn0 + input_temp.in0 + input_temp.jn0 + input_temp.kn0 + input_temp.ln0;
                input_temp.vent_s = input_temp.as0*input_temp.an0 + input_temp.bs0*input_temp.bn0 + input_temp.cs0*input_temp.cn0...
                    + input_temp.ds0*input_temp.dn0 + input_temp.es0*input_temp.en0 + input_temp.fs0*input_temp.fn0+input_temp.gs0;
            end
        case '上下から順次閉じる_面積比から閉じ量求める'
            if (result(i,3) > input.break_pressure || result(i,10) > 0) && result(i,1) <= input.bag_Hs 
                hz = result(i,1);

                if para(12) + para(6) >= result(i,1) && result(i,1) > para(12) - para(6)
                    input_temp.fn0 = input.fn0 * vent_block_ratio3(-hz,-para(12),para(6),'1');

                elseif para(12) - para(6) >= result(i,1)
                    input_temp.fn0 = 0;
                end

                if para(11) + para(5) >= result(i,1) && result(i,1) > para(11) - para(5)
                    input_temp.en0 = input.en0 * vent_block_ratio3(-hz,-para(11),para(5),'1');

                elseif para(11) - para(5) >= result(i,1)
                    input_temp.en0 = 0;
                end

                if para(10) + para(4) >= result(i,1) && result(i,1) > para(10) - para(4)
                    input_temp.dn0 = input.dn0 * vent_block_ratio3(-hz,-para(10),para(4),'1');

                elseif para(10) - para(4) >= result(i,1)
                    input_temp.dn0 = 0;
                end

                if para(9) + para(3) >= result(i,1) && result(i,1) > para(9) - para(3)
                    input_temp.cn0 = input.cn0 * vent_block_ratio3(-hz,-para(9),para(3),'1');

                elseif para(9) - para(3) >= result(i,1)
                    input_temp.cn0 = 0;
                end

                if para(8) + para(2) >= result(i,1) && result(i,1) > para(8) - para(2)
                    input_temp.bn0 = input.bn0 * vent_block_ratio3(-hz,-para(8),para(2),'1');

                elseif para(8) - para(2) >= result(i,1)
                    input_temp.bn0 = 0;
                end

                if para(7) + para(1) >= result(i,1) && result(i,1) > para(7) - para(1)
                    input_temp.an0 = input.an0 * vent_block_ratio3(-hz,-para(7),para(1),'1');

                elseif para(7) - para(1) >= result(i,1)
                    input_temp.an0 = 0;
                end

                input_temp.Vflag  = input_temp.an0 + input_temp.bn0 + input_temp.cn0 + input_temp.dn0 + input_temp.en0 + input_temp.fn0; % + input_temp.gn0 + input_temp.hn0 + input_temp.in0 + input_temp.jn0 + input_temp.kn0 + input_temp.ln0;
                input_temp.vent_s = input_temp.as0*input_temp.an0 + input_temp.bs0*input_temp.bn0 + input_temp.cs0*input_temp.cn0...
                    + input_temp.ds0*input_temp.dn0 + input_temp.es0*input_temp.en0 + input_temp.fs0*input_temp.fn0+input_temp.gs0;
            end
    end

    if Gmax <= result(i,7)
        Gmax = max(result(i,7));
    end

    if Vmin >= abs(result(i,2))
        Vmin = min(abs(result(i,2)));
    end

    % エアバッグがぺったんこになったら強制終了（(；゜Д゜)
    if result(i,1) <= input_temp.bag_dh
        break;
    end

end


end