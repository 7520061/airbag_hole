function state_dot = LandingDynamics(state, input)

global output

% atmos  = input.atmos;
% g = calc1DInterp(atmos.htab, atmos.gtab, state(1,:), atmos.gpp);        % �d�͉����x[m/s2]
% p = calc1DInterp(atmos.htab, atmos.Prestab, state(1,:), atmos.Presspp); % �O�C�� [PaA] (U.S. STANDARD ATMOSPHRE,1976 ���`�⊮)

% �����₷�����邽�߂ɖ��O����
height = state(1);
vel = state(2);
mass = state(3);
tem = state(4);
vol = state(5);

% �����l�̏����ݒ�
d_height = 0;
d_vel = 0;
d_mass = 0;
d_tem = 0;
d_vol = 0;

% �G�A�o�b�O�̏󋵂ɍ��킹�ĉ^�����f����ݒ聕���l��͎��s
if height > input.bag_H  % ���o���h���Ă���Ƃ�
    vertical_falldown_dynamics;
elseif height <= input.bag_H
    vent_landing_dynamics;
elseif  input.bag_p < input.p
    vertical_falldown_dynamics;
end

% �����l�i�[
state_dot = [d_height; d_vel; d_mass; d_tem; d_vol;];

% ���R����
    function vertical_falldown_dynamics
        d_height = vel;
        d_vel = -input.g + (input.d_c(1) * vel + input.d_c(2)) / input.m0;

        % ���ʐ���
        output.bag_h  = input.bag_H;
        output.bag_S  = 0;
        output.bag_p  = input.p;
        output.acc    = -input.g;
        output.G      = (output.acc + input.g) / input.g;
        output.Mae    = 0;
        output.se     = 0;
        output.Te     = 0;
        output.ve     = 0;
        output.rhoe   = 0;
        output.d_mass = 0;
        output.sound_speed = 0;
    end

% �r�C���Ȃ��璅�� ���������̓G�A�o�b�O�̃_�C�i�~�N�X��׋����Ă���i�߂�ׂ��I�I
    function vent_landing_dynamics

%         bag_D = input.bag_D; %�f�ʐϕω��l�����Ȃ�
        bag_D = sqrt(4*input.ex_bag_S/pi); %�����l����f�ʐϕω����l��

        if height>=input.bag_Hs

            r = sqrt((1 - (height - input.bag_Hs)^2 / (input.bag_H0 - input.bag_Hs)^2) * bag_D^2/4);
            bag_S = r^2*pi;
            vol = input.bag_Hs * bag_D^2 * pi / 4 + pi / 3 * (height - input.bag_Hs) * (bag_D^2 / 2 + r^2);

        else
            bag_S = bag_D^2 * pi / 4;
            vol = bag_S * height;

        end
        % �r�C���f���i���k�����̗͊w�̒m�����K�v�@�@�ł��C���Ԃ�]�T�j
        % Mae:�o���}�b�n�� rhoe:�o�����x Pe:�o������
        % ���͂͏�ԕ���������v�Z
        bag_p = mass * input.R * tem / vol;

        if bag_p < input.p
            bag_p = input.p;
        end
        
        vent_se = input.vent_s;
        vent_Mae = sqrt(2 / (input.GAM - 1)*((bag_p / input.p)^((input.GAM - 1) / input.GAM) - 1));
        if vent_Mae <= 0
            vent_Mae = 0;
        elseif 1 < vent_Mae
            vent_Mae = 1;
        end
        if vent_se == 0
            vent_Mae = 0;
        end
        vent_rhoe = mass / vol * (1 + (input.GAM - 1) / 2 * vent_Mae^2)^(-1 / (input.GAM - 1));
        vent_Te = tem * (1 + (input.GAM - 1) / 2 * vent_Mae^2)^(-1);
        sound_speed =  sqrt(input.GAM * input.R * vent_Te);     % ����
        vent_ve = sound_speed * vent_Mae;
        acc = -input.g + (bag_S * (bag_p - input.p) + (input.d_c(1) * vel + input.d_c(2))) / input.m0 ;
        d_height = vel;
        d_vel = acc;
        d_mass =  - vent_rhoe * vent_se * vent_ve  *input.vent_Cd;
        % �������O�C���ƈꏏ�̂Ƃ��̓G�A�o�b�O�����o�E���h���Ă���Ƃ��Ɠ��`�@���Ⴏ�o�b�O�̐ρi�K�X���ʁj�ێ�
        % �����ăo�b�O1�C�����̋K���������X�V
        if bag_p <= input.p && input.bag_V < vol
            input.bag_H = height;
        elseif bag_p <= input.p && input.bag_V >= vol
            d_vol = 0;
            input.bag_H = height;
        else
            d_vol = bag_S * vel;
        end

        if vent_se == 0
            d_tem = -(input.GAM - 1) * tem / vol * d_vol;
        elseif vent_se ~= 0
            d_tem = (input.GAM * input.R * (vent_Te - tem) + 1/2 * vent_ve^2 * (input.GAM - 1)) * d_mass/(mass * input.R * input.GAM);
        end

        %% ���ʐ���
        output.bag_h  = input.bag_H;
        output.bag_S  = bag_S;
        output.bag_p  = bag_p;
        output.acc    = acc;
        output.G      = (acc + input.g) / input.g;
        output.Mae    = vent_Mae;
        output.se     = vent_se;
        output.Te     = vent_Te;
        output.ve     = vent_ve;
        output.rhoe   = vent_rhoe;
        output.d_mass = -d_mass;
        output.sound_speed = sound_speed;
    end
end