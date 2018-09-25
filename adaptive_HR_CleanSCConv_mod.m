function [X, Y, A, mu_store] = adaptive_HR_CleanSCConv_mod(CSM, plane_distance, freqs, ...
    scan_limits, grid_resolution, mic_poses, c, U, Num_sources)

% Beamforming using ADAPTIVE HR-CLEAN-SC according to Luesutthiviboon and
% Sijtsma 2018 adaptive HR-CLEAN-SC paper. Steering vector formulation 
% 'from' Sijtsma. This code INCLUDES convection.
%
% If array is out-of-flow additional SHEAR LAYER CORRECTION has to be
% applied. A simple form is U_corrected = [0 U_wind 0]*(z-z_sl)/z, in this
% case only wind in y-direction, z is beamform distance and z_sl
% approximated shear layer distance from array.
%
% Variables of interest are the:
% - Clean beam width, lambda
% - Loop gain, phi
% - Stop criteria, epsilon
% - Num_sources
%
% Anwar Malgoezar, Sep. 2018
% Group ANCE
%
% Code authors:
% Anwar: Conventional beamforming, CLEAN-SC and convection part
% Salil: HR-CLEAN-SC and adaptive part

    fprintf('\t------------------------------------------\n');
    fprintf('\tStart beamforming, HR CLEAN-SC...\n');

    lambda = 480;
    phi = 0.99; % loop gain
    epsilon = 1e-2;
    % mu = 0.25; % marker location constraint
    max_iter = 20; % maximum iteration number for HR CLEAN SC
    
    % Convection in Mach
    M = U/c;
    beta_sq = 1 - dot(M, M);

    % Setup scanning grid using grid_resolution and dimensions
    N_mic = size(mic_poses, 2);
    N_freqs = length(freqs);

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

    B = zeros(1, N_scanpoints);
    A = zeros(1, N_scanpoints);
    P = zeros(1, N_scanpoints);

    x_0 = mean(mic_poses, 2);
    
    % to store mu
    mu_store = zeros(Num_sources, N_freqs);
    
    r_t0 = sqrt( (x_t(:,1) - x_0(1)).^2 + ...
                 (x_t(:,2) - x_0(2)).^2 + ...
                 (x_t(:,3) - x_0(3)).^2 );
    fprintf('\tBeamforming for %d frequency points...\n', N_freqs);
    for K = 1:N_freqs
        Q = zeros(1, N_scanpoints);

        k = 2*pi*freqs(K)/c;
        h = zeros(N_mic, size(x_t, 1));
        sum_r_ti = zeros(N_scanpoints, 1);
        for I = 1:N_mic
            r_ti_vec = -x_t + mic_poses(:,I).';
            r_ti_sq_adap = sum(M.*r_ti_vec,2).^2 + beta_sq*sum(r_ti_vec.^2,2);
            % Horrible expression
            h(I, :) = exp(-1i*k*( ( (-sum(M.*r_ti_vec,2) + ...
                                     sqrt(r_ti_sq_adap) ) / beta_sq) ...
                                  ) ) ./ (4*pi*sqrt(r_ti_sq_adap));
        end
        
        % Start CLEAN SC procedure
        % Eq. (9), P^(0)

        D = CSM(:,:,K); % To stick to Pieter's conjugate matrix multiplication, I take here the conjugate of the CSM
        CSM_hr = CSM(:,:,K); % Same as D but won't be subtracted by CLEAN-SC
        
        for J = 1:size(h,2)
            g_j = h(:,J);
            P(J) = norm(conj(g_j).'*D.'*g_j)/(norm(g_j))^4;
        end

        % Start iteration
        Cabs = sum(abs(D(:)));
        Dcurr = Cabs;
        count = 0;
        source_location_SC = [];
        Dprev = 1e10;
        while ( Dcurr < Dprev )&&(count<10)
            % Determine peak source
            [Pmax, index_max] = max(abs(P));
            source_location_SC = [source_location_SC, index_max];
            ispositiv = real(P(index_max)) > 0; % Do not 'include' negative pressure maps

            gmax = D.'*(h(:,index_max)./norm(h(:,index_max))^2)/Pmax;
            % Degraded CSM Eq. (14)
            Cource = gmax*conj(gmax).';
            D = D - Pmax*phi*Cource.';
            % Degraded source power Eq. (15), P^(i)
            for J1 = 1:size(h,2)
                g_j = h(:,J1);
                P(J1) = norm(conj(g_j).'*D.'*g_j)/(norm(g_j))^4;
            end

            % Construct clean beam Eq. (13)
            d_frompeak = sqrt( (x_t(:,1) - x_t(index_max,1)).^2 + ...
                               (x_t(:,2) - x_t(index_max,2)).^2 );
            Q = Q + phi*Pmax*10.^(-lambda.*d_frompeak.^2).' * ispositiv;
            
            Dprev = Dcurr;
            Dcurr = sum(abs(D(:)));
            count = count + 1;
        end
        fprintf('\tf = %.2f \n', freqs(K));
        fprintf('\t# of iterations for CLEAN-SC: %d\n', count);

        % Remaining dirty map shouldnt contain negative pressures
        P(real(P)<0) = 0;
        % Add to final map, i.e. sum over frequencies, clean map and residual
        % dirty map
        B = B + Q + P;

        % A = A + P;
        
        %% HR CLEAN-SC
        % obtain number of sources from CLEAN-SC
        if numel(source_location_SC)<Num_sources
            fprintf('\tDetected amount of sources less than set! Num_sources -> %d\n', numel(source_location_SC));
            Num_sources = numel(source_location_SC);
        end
        ind = source_location_SC(1:Num_sources);
        ind_markers_current = ind;
        ind_markers_new = zeros(1,Num_sources);
        ind_sources = zeros(1,Num_sources);
        stop = 0;
        summation_collect = zeros(1,size(h,2));
        count = 0;

        while stop == 0
            A_sum_sources = zeros(1,N_scanpoints);
            CSM_original = CSM_hr;
            mu_iter = zeros(Num_sources, 1);
            for J = 1:Num_sources
                % store mu_check and F_j
                mu_check = zeros(1,size(h,2));
                F_j = zeros(1,size(h,2));
                % source steering vector
                g_j_j = h(:,ind(J));
                for I = 1:size(h,2)
                    g_j = h(:,I);
                    % Equation (7)
                    u_j = g_j./(norm(g_j))^2;
                    mu_check(I) = norm(conj(g_j_j).'*u_j)^2;
                    summation = zeros(size(h,1),1);
                    for L = 1:length(ind)
                        % take into account the other sources only
                        if L ~= J
                            % steering vector for source K
                            g_k = h(:,ind(L));
                            summation = summation + (conj(g_k).'*u_j).*g_k;
                        end
                    end
                    numerator = (norm(summation))^2;
                    summation_collect(:,I) = numerator;
                    % Equation (26)
                    F_j(I) = numerator/(mu_check(I)*(norm(g_j_j))^2);
                end
                
                %%%%%%%%%%% adaptive mu with MSL %%%%%%%%%%%%%%%%%%%
                % plot mu check here
                X_rep = repmat(X,length(Y),1);
                Y_rep = repmat(Y',1,length(X));
                mu_check_plane = reshape(mu_check, N_X, N_Y).';
                %surf(X_rep, Y_rep, mu_check_plane, 'EdgeColor', 'none');
                %view(2)
                %pause(1)
                
                [peaks_x, peaks_y, peaks_mag] = find_peaks(mu_check_plane, X_rep, Y_rep, 0);
                
                if length(peaks_mag) <= 1;
                    mu = min(min(mu_check_plane));
                else
                    mu = peaks_mag(2);
                end
                mu_iter(J,1) = mu;
                
                % find minimum F_j 
                [val_F_j, ind_F_j] = sort(F_j, 'ascend');
                for I = 1:length(ind_F_j)
                    if mu_check(ind_F_j(I)) >= mu
                        ind_minimizer = ind_F_j(I);
                        break
                    end
                end
                ind_markers_new(J) = ind_minimizer;
                % identify the minimizer
                g_j = h(:,ind_minimizer);   % unscaled
                u_j = g_j./(norm(g_j))^2;
                % source component Equation (21)
                % 
                denom = norm(conj(g_j).'*CSM_hr.'*g_j)/(norm(g_j))^4;
                h_j = CSM_hr.'*u_j/denom;
                % estimated source power Equation (25)
                A_j = zeros(1,size(h,2));
                for I = 1:size(h,2)
                    g_i = h(:,I); % unscaled
                    w_i = g_i./(norm(g_i))^2;
                    PSF_i = norm(conj(w_i).'*h_j)^2;
                    A_j(I) = denom*PSF_i;
                end
                % identify the peak in the source map, new source location
                [A_j_max, ind_A_j_max] = max(real(A_j));
                ind_sources(J) = ind_A_j_max;
                d_frompeak = sqrt( (x_t(:,1) - x_t(ind_A_j_max,1)).^2 + ...
                               (x_t(:,2) - x_t(ind_A_j_max,2)).^2 );
                ispositiv = real(A_j_max) > 0; % Do not 'include' negative pressure maps
                
                % degrade CSM
                gmax = D.'*(h(:,ind_A_j_max)./norm(h(:,ind_A_j_max))^2)/A_j_max;
                C_source = gmax*conj(gmax).';
                CSM_hr = CSM_hr - A_j_max*phi*C_source;
                
                % generate a clean beam
                A_sum_sources = A_sum_sources + phi*abs(A_j_max)*10.^(-lambda.*d_frompeak.^2).'*ispositiv;
            end
            % if the markers do not move/ maximum iteration reached, stop
            if isequal(ind_markers_new, ind_markers_current) || count > max_iter
                ind_markers_current = ind_markers_new;
                ind = ind_sources;
%                 mu_store(:,K) = mu_iter;
                % add clean beam plus remaining CSM
                P = sum(conj(h).*(CSM_hr*h), 1);
                P(real(P)<0) = 0;
                A = A + A_sum_sources + P;
                stop = 1;      
            else
                ind_markers_current = ind_markers_new;
                ind = ind_sources;
                CSM_hr = CSM_original;
                count = count + 1;
            end
        end
        fprintf('\t# of iterations for HR CLEAN-SC: %d\n', count);
    end
    % output the estimated source power
    A = reshape(A, N_X, N_Y).';
    fprintf('\tHR CLEAN-SC procedure completed!\n');
    fprintf('\t------------------------------------------\n');
end