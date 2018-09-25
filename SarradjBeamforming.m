function B = SarradjBeamforming(x, CSM, mic_positions, freq, speed_of_sound, steering)
% Beamforming with selection of steering vector. The formulations are based
% on the paper of Sarradj 2012 "Three-Dimensional Acoustic Source Mapping
% with Different Beamforming Steering Vector Formulations". Not the fastest
% code, but good to experiment with.
%
% Anwar Malgoezar, June 2016
% Group ANCE

n_mics = length(mic_positions);
% mic_positions is 3 by n_mic
r = mic_positions - x'*ones(1, n_mics);

R = sqrt( dot(r,r,1) );
xi_timedelay = R/speed_of_sound; 

r_0 = norm(x.' - mean(mic_positions, 2)); % Distance source - array_midpoint
switch steering
    case 'form1'
        h = exp(-1i*2*pi*freq*(R-r_0)/speed_of_sound)/n_mics;
    case 'form2'
        h = R.*exp(-1i*2*pi*freq*(R-r_0)/speed_of_sound)/n_mics/r_0;
    case 'form3'
        h = exp(-1i*2*pi*freq*(R-r_0)/speed_of_sound)./( R*r_0*sum(R.^(-2)) );
    case 'form4'
        h = exp(-1i*2*pi*freq*(R-r_0)/speed_of_sound)./( R.*sqrt(n_mics*sum(R.^(-2))) );
    case 'pieter'
        g = (-exp(-2.*pi.*1i.*freq.*xi_timedelay)) ./ ( 4*pi*R );
        h = g/dot(g,g);
    otherwise
        error('Choose your formulation!');
end

B = h*CSM*h';

end