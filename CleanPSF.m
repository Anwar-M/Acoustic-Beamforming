function [X, Y, B] = CleanPSF(CSM, plane_distance, frequencies, ...
    scan_limits, grid_resolution, mic_positions, c)
% Beamforming using CLEAN PSF according to Sijtsma 2007 CLEAN paper.
% Using steering vector formulation III 'by' Ennes Sarradj. I.e. it 
% is the same as SarradjBeamforming using form3. No diagonal removal.
%
% Anwar Malgoezar, Oct. 2017
% Group ANCE

fprintf('\t------------------------------------------\n');
fprintf('\tStart beamforming, CLEAN-PSF...\n');

lambda = 480;
phi = 0.99; % loop gain

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

fprintf('\tBeamforming for %d frequency points...\n', N_freqs);
for K = 1:N_freqs
    Q = zeros(1, N_scanpoints);
    
    k = 2*pi*frequencies(K)/c;
    h = zeros(N_mic, size(x_t, 1));
    sum_r_ti = zeros(N_scanpoints, 1);
    for I = 1:N_mic
        r_ti = sqrt( (x_t(:,1) - mic_positions(1,I)).^2 + ...
                     (x_t(:,2) - mic_positions(2,I)).^2 + ...
                     (x_t(:,3) - mic_positions(3,I)).^2 );
        sum_r_ti = sum_r_ti + r_ti.^(-2);
        h(I, :) = exp(-1i*k*(r_ti-r_t0))./(r_ti.*r_t0);
    end
    
    for I = 1:N_mic
        h(I, :) = h(I, :) ./ sum_r_ti.';
    end
    
    % Start CLEAN PSF procedure
    % Eq. (9), P^(0)
    P = sum(h.*(CSM(:,:,K)*conj(h)), 1);
    
    % Start iteration
    D = CSM(:,:,K);
    Dprev = 1e10;
    Dcurr = sum(abs(D(:)));
    count = 0;
    while ( Dprev > Dcurr )&&(count < N_mic)
        Dprev = Dcurr;
        
        % Determine peak source
        [Pmax, index_max] = max(abs(P)); % not real() !!
        ispositiv = real(P(index_max)) > 0; % Do not 'include' negative pressure maps
        
        % Steering vector according to maximum point
        gmax = h(:,index_max) * sum_r_ti(index_max) * r_t0(index_max)^2;
        % Degraded CSM Eq. (14)
        D = D - phi*Pmax*conj(gmax)*gmax.';
        % Degraded source power Eq. (15), P^(i)
        P = sum(h.*(D*conj(h)), 1);
        
        % Construct clean beam Eq. (13)
        d_frompeak = sqrt( (x_t(:,1) - x_t(index_max,1)).^2 + ...
                           (x_t(:,2) - x_t(index_max,2)).^2 );
        Q = Q + phi*Pmax*10.^(-lambda.*d_frompeak.^2).' * ispositiv; % Not add negative pressures
        
%         intermediate_check;
        Dcurr = sum(abs(D(:)));
        count = count + 1;
    end
    fprintf('\t# Iterations: %d\n', count);
%     intermediate_check;
    
    % Remaining dirty map shouldnt contain negative pressures
    P(real(P)<0) = 0;
    % Add to final map, i.e. sum over frequencies, clean map and residual
    % dirty map
    B = B + Q + P;
end

B = reshape(B, N_X, N_Y).';

fprintf('\tBeamforming complete!\n');
fprintf('\t------------------------------------------\n');
end