function [t_worst, result_worst, Gmax_worst, Gmax_threesigma, maxload, maxparam] = error_propagation(DynFun, input, t0, state0, para, Gmax)
% シミュレーション頑張るンゴね～

vd_error = 0.0005;
vh_error = 0.001;
bagd_error  = 0.01;
baghs_error = 0.005;

input_temp = input;
state0_temp = state0;

maxload.ad = Gmax; maxload.bd = Gmax; maxload.cd = Gmax;
maxload.dd = Gmax; maxload.ed = Gmax; maxload.fd = Gmax;

maxload.ah = Gmax; maxload.bh = Gmax; maxload.ch = Gmax;
maxload.dh = Gmax; maxload.eh = Gmax; maxload.fh = Gmax;

maxload.baghs = Gmax; maxload.bagd = Gmax;

maxparam.ad = para(1);  maxparam.bd = para(2);  maxparam.cd = para(3);
maxparam.dd = para(4);  maxparam.ed = para(5);  maxparam.fd = para(6);

maxparam.ah = para(7);  maxparam.bh = para(8);  maxparam.ch = para(9);
maxparam.dh = para(10); maxparam.eh = para(11); maxparam.fh = para(12);

maxparam.bagd = input.bag_D;  maxparam.baghs = input.bag_Hs;

for vd_err = 1 : 1 : 6
    paratemp  = para;
    if para(vd_err) ~= 0.0

        paratemp(vd_err) = para(vd_err) + vd_error;

        input_temp.as0 = paratemp(1)^2 / 4 * pi * input.an0;
        input_temp.bs0 = paratemp(2)^2 / 4 * pi * input.bn0;
        input_temp.cs0 = paratemp(3)^2 / 4 * pi * input.cn0;
        input_temp.ds0 = paratemp(4)^2 / 4 * pi * input.dn0;
        input_temp.es0 = paratemp(5)^2 / 4 * pi * input.en0;
        input_temp.fs0 = paratemp(6)^2 / 4 * pi * input.fn0;

        [~, ~, Gmax_temp] = AirbagSim(DynFun, input_temp, t0, state0, paratemp);

        if vd_err == 1
            if Gmax_temp > Gmax
                maxload.ad = Gmax_temp;
                maxparam.ad = paratemp(vd_err);
            end
        elseif vd_err == 2
            if Gmax_temp > Gmax
                maxload.bd = Gmax_temp;
                maxparam.bd = paratemp(vd_err);
            end
        elseif vd_err == 3
            if Gmax_temp > Gmax
                maxload.cd = Gmax_temp;
                maxparam.cd = paratemp(vd_err);
            end
        elseif vd_err == 4
            if Gmax_temp > Gmax
                maxload.dd = Gmax_temp;
                maxparam.dd = paratemp(vd_err);
            end
        elseif vd_err == 5
            if Gmax_temp > Gmax
                maxload.ed = Gmax_temp;
                maxparam.ed = paratemp(vd_err);
            end
        elseif vd_err == 6
            if Gmax_temp > Gmax
                maxload.fd = Gmax_temp;
                maxparam.fd = paratemp(vd_err);
            end
        end
        paratemp  = para;
        paratemp(vd_err) = para(vd_err) - vd_error;

        input_temp.as0 = paratemp(1)^2 / 4 * pi * input.an0;
        input_temp.bs0 = paratemp(2)^2 / 4 * pi * input.bn0;
        input_temp.cs0 = paratemp(3)^2 / 4 * pi * input.cn0;
        input_temp.ds0 = paratemp(4)^2 / 4 * pi * input.dn0;
        input_temp.es0 = paratemp(5)^2 / 4 * pi * input.en0;
        input_temp.fs0 = paratemp(6)^2 / 4 * pi * input.fn0;

        [~, ~, Gmax_temp] = AirbagSim(DynFun, input_temp, t0, state0, paratemp);

        if vd_err == 1
            if Gmax_temp > Gmax
                maxload.ad = Gmax_temp;
                maxparam.ad = paratemp(vd_err);
            end
        elseif vd_err == 2
            if Gmax_temp > Gmax
                maxload.bd = Gmax_temp;
                maxparam.bd = paratemp(vd_err);
            end
        elseif vd_err == 3
            if Gmax_temp > Gmax
                maxload.cd = Gmax_temp;
                maxparam.cd = paratemp(vd_err);
            end
        elseif vd_err == 4
            if Gmax_temp > Gmax
                maxload.dd = Gmax_temp;
                maxparam.dd = paratemp(vd_err);
            end
        elseif vd_err == 5
            if Gmax_temp > Gmax
                maxload.ed = Gmax_temp;
                maxparam.ed = paratemp(vd_err);
            end
        elseif vd_err == 6
            if Gmax_temp > Gmax
                maxload.fd = Gmax_temp;
                maxparam.fd = paratemp(vd_err);
            end
        end
    end
end

input_temp.fs0 = para(6)^2 / 4 * pi * input.fn0;

for vh_err = 7 : 1 : 12
    paratemp  = para;
    if para(vh_err-6) ~= 0.0
        paratemp(vh_err) = para(vh_err) + vh_error;

        [~, ~, Gmax_temp] = AirbagSim(DynFun, input_temp, t0, state0, paratemp);

        if vh_err == 7
            if Gmax_temp > Gmax
                maxload.ah = Gmax_temp;
                maxparam.ah = paratemp(vh_err);
            end
        elseif vh_err == 8
            if Gmax_temp > Gmax
                maxload.bh = Gmax_temp;
                maxparam.bh = paratemp(vh_err);
            end
        elseif vh_err == 9
            if Gmax_temp > Gmax
                maxload.ch = Gmax_temp;
                maxparam.ch = paratemp(vh_err);
            end
        elseif vh_err == 10
            if Gmax_temp > Gmax
                maxload.dh = Gmax_temp;
                maxparam.dh = paratemp(vh_err);
            end
        elseif vh_err == 11
            if Gmax_temp > Gmax
                maxload.eh = Gmax_temp;
                maxparam.eh = paratemp(vh_err);
            end
        elseif vh_err == 12
            if Gmax_temp > Gmax
                maxload.fh = Gmax_temp;
                maxparam.fh = paratemp(vh_err);
            end
        end
        paratemp  = para;
        paratemp(vh_err) = para(vh_err) - vh_error;

        [~, ~, Gmax_temp] = AirbagSim(DynFun, input_temp, t0, state0, paratemp);

        if vh_err == 7
            if Gmax_temp > Gmax
                maxload.ah = Gmax_temp;
                maxparam.ah = paratemp(vh_err);
            end
        elseif vh_err == 8
            if Gmax_temp > Gmax
                maxload.bh = Gmax_temp;
                maxparam.bh = paratemp(vh_err);
            end
        elseif vh_err == 9
            if Gmax_temp > Gmax
                maxload.ch = Gmax_temp;
                maxparam.ch = paratemp(vh_err);
            end
        elseif vh_err == 10
            if Gmax_temp > Gmax
                maxload.dh = Gmax_temp;
                maxparam.dh = paratemp(vh_err);
            end
        elseif vh_err == 11
            if Gmax_temp > Gmax
                maxload.eh = Gmax_temp;
                maxparam.eh = paratemp(vh_err);
            end
        elseif vh_err == 12
            if Gmax_temp > Gmax
                maxload.fh = Gmax_temp;
                maxparam.fh = paratemp(vh_err);
            end
        end
    end
end

paratemp  = para;

input_temp.bag_D = input.bag_D + bagd_error;

input_temp.bag_V = 4 / 3 * pi * input_temp.bag_D^2 / 4 * (input.bag_H0 - input.bag_Hs) / 2 + input_temp.bag_D^2 * pi / 4 * input.bag_Hs; % バッグ体積(ガス量過不足なし)[m3]
input_temp.mol    = input.p0 * input_temp.bag_V /(input.R0 * input.T0); % 充填ガスモル量[mol]
input_temp.gas_m0 = (input_temp.mol * input.atomicm - input.gas_margin) * 10^-3;

input_temp.bag_V0 = input_temp.mol * input.R0 * input.T0 / input.p0; % 初期バッグ体積　[m3], ゲージ圧0, ガス量に合わせて体積は小さくなる計算
input_temp.bag_rho0 = input_temp.gas_m0 / input_temp.bag_V0;  % 初期密度[kg/m3]

state0_temp(3) = input_temp.gas_m0;
state0_temp(5) = input_temp.bag_V0;

[~, ~, Gmax_temp] = AirbagSim(DynFun, input_temp, t0, state0_temp, paratemp);
if Gmax_temp > Gmax
    maxload.bagd = Gmax_temp;
    maxparam.bagd = input_temp.bag_D;
end

input_temp.bag_D = input.bag_D - bagd_error;

input_temp.bag_V = 4 / 3 * pi * input_temp.bag_D^2 / 4 * (input.bag_H0 - input.bag_Hs) / 2 + input_temp.bag_D^2 * pi / 4 * input.bag_Hs; % バッグ体積(ガス量過不足なし)[m3]
input_temp.mol    = input.p0 * input_temp.bag_V /(input.R0 * input.T0); % 充填ガスモル量[mol]
input_temp.gas_m0 = (input_temp.mol * input.atomicm - input.gas_margin) * 10^-3;

input_temp.bag_V0 = input_temp.mol * input.R0 * input.T0 / input.p0; % 初期バッグ体積　[m3], ゲージ圧0, ガス量に合わせて体積は小さくなる計算
input_temp.bag_rho0 = input_temp.gas_m0 / input_temp.bag_V0;  % 初期密度[kg/m3]

state0_temp(3) = input_temp.gas_m0;
state0_temp(5) = input_temp.bag_V0;

[~, ~, Gmax_temp] = AirbagSim(DynFun, input_temp, t0, state0_temp, paratemp);
if Gmax_temp > maxload.bagd
    maxload.bagd = Gmax_temp;
    maxparam.bagd = input_temp.bag_D;
end

input_temp.bag_D = input.bag_D;

input_temp.bag_Hs = input.bag_Hs + baghs_error;
input_temp.bag_H0 = 0.25 + input_temp.bag_Hs;  % 展開時バッグ高さ [m]

input_temp.h0 = input_temp.bag_H0 + 0.001;% + 2; % Initial Altitude @ground level [m]

input_temp.bag_V = 4 / 3 * pi * input_temp.bag_D^2 / 4 * (input_temp.bag_H0 - input_temp.bag_Hs) / 2 + input_temp.bag_D^2 * pi / 4 * input_temp.bag_Hs; % バッグ体積(ガス量過不足なし)[m3]
input_temp.mol    = input.p0 * input_temp.bag_V /(input.R0 * input.T0); % 充填ガスモル量[mol]
input_temp.gas_m0 = (input_temp.mol * input.atomicm - input.gas_margin) * 10^-3;

input_temp.bag_V0 = input_temp.mol * input.R0 * input.T0 / input.p0; % 初期バッグ体積　[m3], ゲージ圧0, ガス量に合わせて体積は小さくなる計算
input_temp.bag_rho0 = input_temp.gas_m0 / input_temp.bag_V0;  % 初期密度[kg/m3]

state0_temp(1) = input_temp.h0;
state0_temp(3) = input_temp.gas_m0;
state0_temp(5) = input_temp.bag_V0;

[~, ~, Gmax_temp] = AirbagSim(DynFun, input_temp, t0, state0_temp, paratemp);
if Gmax_temp > Gmax
    maxload.baghs = Gmax_temp;
    maxparam.baghs = input_temp.bag_Hs;
end

input_temp.bag_Hs = input.bag_Hs - baghs_error;
input_temp.bag_H0 = 0.25 + input_temp.bag_Hs;  % 展開時バッグ高さ [m]

input_temp.h0 = input_temp.bag_H0 + 0.001;% + 2; % Initial Altitude @ground level [m]

input_temp.bag_V = 4 / 3 * pi * input_temp.bag_D^2 / 4 * (input_temp.bag_H0 - input_temp.bag_Hs) / 2 + input_temp.bag_D^2 * pi / 4 * input_temp.bag_Hs; % バッグ体積(ガス量過不足なし)[m3]
input_temp.mol    = input.p0 * input_temp.bag_V /(input.R0 * input.T0); % 充填ガスモル量[mol]
input_temp.gas_m0 = (input_temp.mol * input.atomicm - input.gas_margin) * 10^-3;

input_temp.bag_V0 = input_temp.mol * input.R0 * input.T0 / input.p0; % 初期バッグ体積　[m3], ゲージ圧0, ガス量に合わせて体積は小さくなる計算
input_temp.bag_rho0 = input_temp.gas_m0 / input_temp.bag_V0;  % 初期密度[kg/m3]

state0_temp(1) = input_temp.h0;
state0_temp(3) = input_temp.gas_m0;
state0_temp(5) = input_temp.bag_V0;

[~, ~, Gmax_temp] = AirbagSim(DynFun, input_temp, t0, state0_temp, paratemp);
if Gmax_temp > maxload.baghs
    maxload.baghs = Gmax_temp;
    maxparam.baghs = input_temp.bag_Hs;
end

paratemp(1) = maxparam.ad; paratemp(2) = maxparam.bd; paratemp(3) = maxparam.cd;
paratemp(4) = maxparam.dd; paratemp(5) = maxparam.ed; paratemp(6) = maxparam.fd;
%         input_temp.gd0 = maxparam.ad; input_temp.hd0 = maxparam.bd; input_temp.id0 = maxparam.cd;
%         input_temp.jd0 = maxparam.dd; input_temp.kd0 = maxparam.ed; input_temp.ld0 = maxparam.fd;

input_temp.as0 = paratemp(1)^2 / 4 * pi * input.an0;
input_temp.bs0 = paratemp(2)^2 / 4 * pi * input.bn0;
input_temp.cs0 = paratemp(3)^2 / 4 * pi * input.cn0;
input_temp.ds0 = paratemp(4)^2 / 4 * pi * input.dn0;
input_temp.es0 = paratemp(5)^2 / 4 * pi * input.en0;
input_temp.fs0 = paratemp(6)^2 / 4 * pi * input.fn0;

paratemp(7)  = maxparam.ah; paratemp(8)  = maxparam.bh; paratemp(9)  = maxparam.ch;
paratemp(10) = maxparam.dh; paratemp(11) = maxparam.eh; paratemp(12) = maxparam.fh;
%         input_temp.gh0 = maxparam.ah; input_temp.hh0 = maxparam.bh; input_temp.ih0 = maxparam.ch;
%         input_temp.jh0 = maxparam.dh; input_temp.kh0 = maxparam.eh; input_temp.lh0 = maxparam.fh;

input_temp.bag_D = maxparam.bagd;
input_temp.bag_Hs = maxparam.baghs; input_temp.bag_H0 = 0.25 + input_temp.bag_Hs;

input_temp.h0 = input_temp.bag_H0 + 0.001;% + 2; % Initial Altitude @ground level [m]

input_temp.bag_V = 4 / 3 * pi * input_temp.bag_D^2 / 4 * (input_temp.bag_H0 - input_temp.bag_Hs) / 2 + input_temp.bag_D^2 * pi / 4 * input_temp.bag_Hs; % バッグ体積(ガス量過不足なし)[m3]
input_temp.mol    = input.p0 * input_temp.bag_V /(input.R0 * input.T0); % 充填ガスモル量[mol]
input_temp.gas_m0 = (input_temp.mol * input.atomicm - input.gas_margin) * 10^-3;

input_temp.bag_V0 = input_temp.mol * input.R0 * input.T0 / input.p0; % 初期バッグ体積　[m3], ゲージ圧0, ガス量に合わせて体積は小さくなる計算
input_temp.bag_rho0 = input_temp.gas_m0 / input_temp.bag_V0;  % 初期密度[kg/m3]

state0_temp(1) = input_temp.h0;
state0_temp(3) = input_temp.gas_m0;
state0_temp(5) = input_temp.bag_V0;

[t_worst, result_worst, Gmax_temp] = AirbagSim(DynFun, input_temp, t0, state0_temp, paratemp);

Gmax_worst = Gmax_temp;

if maxparam.ad - para(1) ~= 0
    dfdx_ad = (maxload.ad - Gmax)/(maxparam.ad - para(1));
else
    dfdx_ad = 0;
end
if maxparam.bd - para(2) ~= 0
    dfdx_bd = (maxload.bd - Gmax)/(maxparam.bd - para(2));
else
    dfdx_bd = 0;
end
if maxparam.cd - para(3) ~= 0
    dfdx_cd = (maxload.cd - Gmax)/(maxparam.cd - para(3));
else
    dfdx_cd = 0;
end
if maxparam.dd - para(4) ~= 0
    dfdx_dd = (maxload.dd - Gmax)/(maxparam.dd - para(4));
else
    dfdx_dd = 0;
end
if maxparam.ed - para(5) ~= 0
    dfdx_ed = (maxload.ed - Gmax)/(maxparam.ed - para(5));
else
    dfdx_ed = 0;
end
if maxparam.fd - para(6) ~= 0
    dfdx_fd = (maxload.fd - Gmax)/(maxparam.fd - para(6));
else
    dfdx_fd = 0;
end

if maxparam.ah - para(7) ~= 0
    dfdx_ah = (maxload.ah - Gmax)/(maxparam.ah - para(7));
else
    dfdx_ah = 0;
end
if maxparam.bh - para(8) ~= 0
    dfdx_bh = (maxload.bh - Gmax)/(maxparam.bh - para(8));
else
    dfdx_bh = 0;
end
if maxparam.ch - para(9) ~= 0
    dfdx_ch = (maxload.ch - Gmax)/(maxparam.ch - para(9));
else
    dfdx_ch = 0;
end
if maxparam.dh - para(10) ~= 0
    dfdx_dh = (maxload.dh - Gmax)/(maxparam.dh - para(10));
else
    dfdx_dh = 0;
end
if maxparam.eh - para(11) ~= 0
    dfdx_eh = (maxload.eh - Gmax)/(maxparam.eh - para(11));
else
    dfdx_eh = 0;
end
if maxparam.fh - para(12) ~= 0
    dfdx_fh = (maxload.fh - Gmax)/(maxparam.fh - para(12));
else
    dfdx_fh = 0;
end

if maxparam.bagd - input.bag_D ~= 0
    dfdx_bagd = (maxload.bagd - Gmax)/(maxparam.bagd - input.bag_D);
else
    dfdx_bagd = 0;
end
if maxparam.baghs - input.bag_Hs ~= 0
    dfdx_baghs = (maxload.baghs - Gmax)/(maxparam.baghs - input.bag_Hs);
else
    dfdx_baghs = 0;
end

sigmax.ad = (maxparam.ad - para(1))/3;
sigmax.bd = (maxparam.bd - para(2))/3;
sigmax.cd = (maxparam.cd - para(3))/3;
sigmax.dd = (maxparam.dd - para(4))/3;
sigmax.ed = (maxparam.ed - para(5))/3;
sigmax.fd = (maxparam.fd - para(6))/3;
sigmax.ah = (maxparam.ah - para(7))/3;
sigmax.bh = (maxparam.bh - para(8))/3;
sigmax.ch = (maxparam.ch - para(9))/3;
sigmax.dh = (maxparam.dh - para(10))/3;
sigmax.eh = (maxparam.eh - para(11))/3;
sigmax.fh = (maxparam.fh - para(12))/3;
sigmax.bagd = (maxparam.bagd - input.bag_D)/3;
sigmax.baghs = (maxparam.baghs - input.bag_Hs)/3;

sigma_y = sqrt((dfdx_ad^2 * sigmax.ad^2 + dfdx_bd^2 * sigmax.bd^2 + dfdx_cd^2 * sigmax.cd^2 + dfdx_dd^2 * sigmax.dd^2 + dfdx_ed^2 * sigmax.ed^2 + dfdx_fd^2 * sigmax.fd^2 + dfdx_ah^2 * sigmax.ah^2 + dfdx_bh^2 * sigmax.bh^2 + dfdx_ch^2 * sigmax.ch^2 + dfdx_dh^2 * sigmax.dh^2 + dfdx_eh^2 * sigmax.eh^2 + dfdx_fh^2 * sigmax.fh^2) * 2 + dfdx_bagd^2 * sigmax.bagd^2 + dfdx_baghs^2 * sigmax.baghs^2);
Gmax_threesigma = Gmax + 3 * sigma_y;

%     fprintf('\nGmax = %f\n',Gmax_temp)


end



