function [peaks_x, peaks_y, peaks_mag] = find_peaks(surface, scan_plane_x, scan_plane_y, rank)
    peaks_x = []; peaks_y = []; peaks_mag = [];
    [m, n] = size(surface);
    for i = 2:m-1;
        for j = 2:n-1;
            test_peak = surface(i,j);
            test_matrix = [1,0;1,1;0,1;-1,1;-1,0;-1,-1;0,-1;1,-1];
            k = 1;
            while k <= 8;
                test_neighbor = surface(i+test_matrix(k,1),j+test_matrix(k,2));
                if test_neighbor > test_peak;
                    break
                end
                if k == 8;
                    [real_x, real_y] = get_real_coordinates_original([i; j], scan_plane_x, scan_plane_y); 
                    peaks_x = [peaks_x; real_x];
                    peaks_y = [peaks_y; real_y];
                    peaks_mag = [peaks_mag; test_peak];
                end
                k = k + 1;
            end
        end
    end
    [peaks_mag, peaks_order] = sort(peaks_mag,'descend');
    peaks_x = peaks_x(peaks_order,:);
    peaks_y = peaks_y(peaks_order,:);
    if rank ~= 0;
        peaks_mag = peaks_mag(1:rank);
        peaks_x = peaks_x(1:rank);
        peaks_y = peaks_y(1:rank);
    end
    
    
    
            
    