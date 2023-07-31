function [atmos] = loadAtmos_airbag()

% U.S. 1976 Standard Atmosphere Table
% Reference:  U.S. 1976 Standard Atmosphere, National Oceanographic
% and Atmospheric Administration, 1976.
% Column 1:  Geometric Altitude (m)
% Column 2:  Atmospheric Density (kg/m^3) Table I
% Column 3:  Speed of Sound (m/s) Table III
% Column 4:  Temperature (K) Table I
% Column 5:  Pressure (Pa) Table I
% Column 6:  Gravity (m/s2) Table II
us1976 = [-2000  1.4782e+00  3.4789e+02  301.154   1.2778e+05  9.8128
              0  1.2250e+00  3.4029e+02  288.150  1.01325e+05  9.8066
             50  1.2191e+00  3.4010e+02  287.825   1.0072e+05  9.8065
            300  1.1901e+00  3.3914e+02  286.200   9.7772e+04  9.8057
            350  1.1844e+00  3.3895e+02  285.875   9.7190e+04  9.8056
           1000  1.1117e+00  3.3643e+02  281.651   8.9876e+04  9.8036
           1500  1.0581e+00  3.3449e+02  278.402   8.4559e+04  9.8020
           2000  1.0066e+00  3.3253e+02  275.154   7.9501e+04  9.8005
           4000  8.1935e-01  3.2459e+02  262.166   6.1660e+04  9.7943];
        
GAM = 1.402;
R = 287.04;
us1976(:,3) = sqrt(GAM * R * us1976(:,4));  % Speed of sound is calculated assuming the constant gamma

atmos.Prestab = us1976(:,5); % us1976(:,2).*287.*us1976(:,4);  % From state equations of air  
atmos.htab    = us1976(:,1);
atmos.Rhotab  = us1976(:,2);
atmos.SOStab  = us1976(:,3);
atmos.gtab    = us1976(:,6);

atmos.Rhopp   = griddedInterpolant(atmos.htab, atmos.Rhotab, 'linear');
atmos.SOSpp   = griddedInterpolant(atmos.htab, atmos.SOStab, 'linear');
atmos.Presspp = griddedInterpolant(atmos.htab, atmos.Prestab, 'linear');
atmos.gpp     = griddedInterpolant(atmos.htab, atmos.gtab, 'linear');

% atmos.Rhopp = interp1(us1976(:,1),us1976(:,2),'pchip','pp');
% atmos.SOSpp = interp1(us1976(:,1),us1976(:,3),'pchip','pp');
% % atmos.Temppp = interp1(us1976(:,1),us1976(:,4),'pchip','pp');
% % atmos.Prespp = interp1(us1976(:,1), us1976(:,2).*287.*us1976(:,4), 'pchip', 'pp');  % From state equations of air  

end