function [X, Y, B] = FastBeamforming3Conv(CSM, plane_distance, frequencies, ...
    scan_limits, grid_resolution, mic_positions, c, U)
% Fast (fewer loops, using array manipulation) beamforming considering
% steering vector formulation III 'by' Ennes Sarradj. I.e. it is the same
% as SarradjBeamforming using form3. This code is faster, but less
% intuitive to understand. This code INCLUDES convection.
%
% If array is out-of-flow additional SHEAR LAYER CORRECTION has to be
% applied. A simple form is U_corrected = [0 U_wind 0]*(z-z_sl)/z, in this
% case only wind in y-direction, z is beamform distance and z_sl
% approximated shear layer distance from array.

% fprintf('\t------------------------------------------\n');
%
% Anwar Malgoezar, Nov. 2017
% Group ANCE

fprintf('\tStart beamforming...\n');

% Convection in Mach
M = U/c;
beta_sq = 1 - dot(M, M);

% Setup scanning grid using grid_resolution and dimensions
N_mic = size(mic_positions, 2);
N_freqs = length(frequencies);

X = scan_limits(1):grid_resolution:scan_limits(2);
Y = scan_limits(3):grid_resolution:scan_limits(4);
Z = plane_distance;
N_X = length(X);
N_Y = length(Y);
N_Z = length(Z);

N_scanpoints = N_X*N_Y*N_Z;
x_t = zeros(N_scanpoints, 3);

x_t(:, 1) = repmat(X, 1, N_Y);
dummy = repmat(Y, N_X, 1);
x_t(:, 2) = dummy(:);
x_t(:, 3) = plane_distance;

x_0 = mean(mic_positions, 2);
B = zeros(1, N_scanpoints);

r_t0 = sqrt( (x_t(:,1) - x_0(1)).^2 + ...
             (x_t(:,2) - x_0(2)).^2 + ...
             (x_t(:,3) - x_0(3)).^2 );

reverseStr = '';
for K = 1:N_freqs
    msg = sprintf('\tBeamforming %d/%d frequency points...\n', K, N_freqs);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    k = 2*pi*frequencies(K)/c;
    h = zeros(N_mic, size(x_t, 1));
    sum_r_ti = zeros(N_scanpoints, 1);
    for I = 1:N_mic
        r_ti_vec = -x_t + mic_positions(:,I).';
        r_ti_sq_adap = sum(M.*r_ti_vec,2).^2 + beta_sq*sum(r_ti_vec.^2,2);
        sum_r_ti = sum_r_ti + 1./r_ti_sq_adap;
        % Horrible expression, evaluation at ARRAY CENTER
        h(I, :) = exp(-1i*k*( ( (-sum(M.*r_ti_vec,2) + ...
                                 sqrt(r_ti_sq_adap) ) / beta_sq) ...
                             - r_t0 ) ) ./ (sqrt(r_ti_sq_adap).*r_t0);
    end
    
    for I = 1:N_mic
        h(I, :) = h(I, :) ./ sum_r_ti.';
    end
    
    B = B + sum(h.*(CSM(:,:,K)*conj(h)), 1);

end

B = reshape(B, N_X, N_Y).';

fprintf('\tBeamforming complete!\n');
fprintf('\t------------------------------------------\n');
end