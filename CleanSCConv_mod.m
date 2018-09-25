function [X, Y, B] = CleanSCConv_mod(CSM, plane_distance, frequencies, ...
    scan_limits, grid_resolution, mic_positions, c, U)
% WEIRD BEHAVIOUR! Beamforming using CLEAN SC according to Sijtsma 2007 CLEAN paper.
% Using steering vector of Pieter Sijtsma. No diagonal removal.
% Variables of interest are the:
% - Clean beam width, lambda
% - Loop gain, phi
% - Stop criteria, epsilon
%
% Anwar Malgoezar, Aug. 2018
% Group ANCE

fprintf('\t------------------------------------------\n');
fprintf('\tStart beamforming, CLEAN-SC...\n');

lambda = 480;
phi = 0.99; % loop gain
epsilon = 1e-2;

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

fprintf('\tBeamforming for %d frequency points...\n', N_freqs);
for K = 1:N_freqs
    Q = zeros(1, N_scanpoints);
    
    k = 2*pi*frequencies(K)/c;
    h = zeros(N_mic, size(x_t, 1));
    sum_r_ti = zeros(N_scanpoints, 1);
    for I = 1:N_mic
        r_ti_vec = -x_t + mic_positions(:,I).';
        r_ti_sq_adap = sum(M.*r_ti_vec,2).^2 + beta_sq*sum(r_ti_vec.^2,2);
        % Horrible expression
        sum_r_ti = sum_r_ti + 1./r_ti_sq_adap;
        h(I, :) = 4*pi*exp(-1i*k*( ( (-sum(M.*r_ti_vec,2) + ...
                                 sqrt(r_ti_sq_adap) ) / beta_sq) ...
                                 ) ) ./ (sqrt(r_ti_sq_adap));
    end
    
    for I = 1:N_mic
        h(I, :) = h(I, :) ./ sum_r_ti.';
    end
    
    % Start CLEAN SC procedure
    % Eq. (9), P^(0)
    
    D = conj(CSM(:,:,K)); % To stick to Pieter's conjugate matrix multiplication, I take here the conjugate of the CSM
    P = sum(conj(h).*(D*h), 1);
    
    % Start iteration
    Cabs = sum(abs(D(:)));
    Dcurr = Cabs;
    count = 0;
    Dprev = 1e10;
    while ( Dcurr < Dprev )&&(count<10)
%     while ( Dcurr > epsilon*Cabs )
        % Determine peak source
        [Pmax, index_max] = max(abs(P));
        ispositiv = real(P(index_max)) > 0; % Do not 'include' negative pressure maps
        
        gmax = D*h(:,index_max)/Pmax;
        % Degraded CSM Eq. (14)
        Cource = Pmax*(gmax)*gmax';
        D = D - phi*Cource;
        % Degraded source power Eq. (15), P^(i)
        P = sum(conj(h).*(D*h), 1);
        
        % Construct clean beam Eq. (13)
        d_frompeak = sqrt( (x_t(:,1) - x_t(index_max,1)).^2 + ...
                           (x_t(:,2) - x_t(index_max,2)).^2 );
        Q = Q + phi*Pmax*10.^(-lambda.*d_frompeak.^2).' * ispositiv;
        
        Dprev = Dcurr;
        Dcurr = sum(abs(D(:)));
        count = count + 1;
%         intermediate_check;
    end
    fprintf('\tf = %.2f, # of iterations: %d\n', frequencies(K), count);
    
    % Remaining dirty map shouldnt contain negative pressures
    P(real(P)<0) = 0;
    % Add to final map, i.e. sum over frequencies, clean map and residual
    % dirty map
    B = B + Q + P;
end

B = reshape(B, N_X, N_Y).';

fprintf('\tCLEAN-SC procedure completed!\n');
fprintf('\t------------------------------------------\n');
end