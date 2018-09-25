function [real_x, real_y] = get_real_coordinates_original(vector_ij, scan_plane_x, scan_plane_y)
    x_1_ind = floor(vector_ij(1,1));
    x_2_ind = ceil(vector_ij(1,1));
    x_offset_norm = vector_ij(1,1) - x_1_ind;
    y_1_ind = floor(vector_ij(2,1));
    y_2_ind = ceil(vector_ij(2,1));
    y_offset_norm = vector_ij(2,1) - y_1_ind;
    % check something here?
    x_1 = scan_plane_x(x_1_ind, y_1_ind);
    x_2 = scan_plane_x(x_2_ind, y_1_ind);
    y_1 = scan_plane_y(x_1_ind, y_1_ind);
    y_2 = scan_plane_y(x_1_ind, y_2_ind);
    x_diff = x_2 - x_1;
    y_diff = y_2 - y_1;
    real_x = x_1 + (x_offset_norm*x_diff);
    real_y = y_1 + (y_offset_norm*y_diff);