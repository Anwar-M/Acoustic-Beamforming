function [X, Y, B] = FastBeamformingFunc(CSM, plane_distance, frequencies, ...
    scan_limits, grid_resolution, mic_positions, c, nu)
% Fast (fewer loops, using array manipulation) beamforming considering
% steering vector formulation III 'by' Ennes Sarradj. I.e. it is the same
% as SarradjBeamforming.m using form3. This code is faster, but less
% intuitive to understand. Functional, not really working well.

% fprintf('\t------------------------------------------\n');
%
% Anwar Malgoezar, June 2018
% Group ANCE

fprintf('\tStart beamforming...\n');

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
    
%     CSM(:,:,K)=reshape(CSM(:,:,K),N_mic,N_mic); 
    CSMopt = balance(CSM(:,:,K),'noperm');
    [U,sigma] = eig(CSMopt);
    sigma(sigma<0)=0;
    CSM_new = U*(sigma.^(1/nu))*U';
    
    msg = sprintf('\tBeamforming %d/%d frequency points...\n', K, N_freqs);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    k = 2*pi*frequencies(K)/c;
    h = zeros(N_mic, size(x_t, 1));
    sum_r_ti = zeros(N_scanpoints, 1);
    for I = 1:N_mic
        r_ti = sqrt( (x_t(:,1) - mic_positions(1,I)).^2 + ...
                     (x_t(:,2) - mic_positions(2,I)).^2 + ...
                     (x_t(:,3) - mic_positions(3,I)).^2 );
        sum_r_ti = sum_r_ti + r_ti.^(-2);
        h(I, :) = -exp(-1i*k*r_ti)./(4*pi*r_ti);
    end
    h = h./dot(h,h,1);

    B = B + sum(h.*(CSM_new*conj(h)), 1).^(nu);
%     B = B + sum(((h.*(CSM_new(:,:,K)*conj(h))).^nu)./((sum_r_ti.').^(2*nu)), 1);

end

B = reshape(B, N_X, N_Y).';

fprintf('\tBeamforming complete!\n');
fprintf('\t------------------------------------------\n');
end