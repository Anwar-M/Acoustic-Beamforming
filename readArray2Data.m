function [DATA, ARRAY_MICS] = readArray2Data(file_path, remove_mic, calib, calib_file_path)
%   ---------------------------------------------------
%   Reading array data for TU Delft acoustic arrays
%   Author: Anwar Malgoezar
%   Date: 2-12-2016
%	---------------------------------------------------

fprintf('\t------------------------------------------\n');
fprintf('\tRead in measurement data from TU Delft acoustic array 2\n');

% LabVIEW calibrated data or MATLAB
if calib
    fprintf('\tLoading LabVIEW calibrated data...\n');
    fid = fopen([file_path '\acoustic_data.cal'],'r');
    DATA = fread(fid, [64 inf], 'single', 0, 'l').';
    fclose(fid);
    fprintf('\tLabVIEW calibrated data loaded!\n');
else
    fprintf('\tLoading RAW data...\n');
    fid = fopen([file_path '\acoustic_data'],'r');
    DATA = fread(fid,[64 inf],'int16',0,'l').';
    fclose(fid);
    fprintf('\tPerforming simple MATLAB calibration...\n');
    
    % Get amplification factor from info.txt
    fid = fopen([file_path '\info.txt'],'r');
    fgetl(fid);fgetl(fid);fgetl(fid);fgetl(fid);fgetl(fid);
    amp_line = fgetl(fid);
    fclose(fid);
    daq_amp_factor = (amp_line(15)=='L')*3.985+(amp_line(15)=='H')*28.03;
    
    % Perform simple calibration
    DATA = DATA*2.5/32768/daq_amp_factor;
    mic_responses = importdata(calib_file_path,'\t',0);
    mic_responses(:, [1 3 4]) = [];
    
    array2mic_order = 8*(1:64) + 56 - 63 * ceil((1:64)/8);
    DATA(:, array2mic_order) = DATA(:, array2mic_order)./(ones(size(DATA,1),1)*mic_responses(1:64).');
    DATA = DATA - (ones(size(DATA,1),1)*mean(DATA));
    fprintf('\tRaw data calibrated!\n');
end

% Load microphone positions
fid = fopen([file_path '\configuration\config.txt'],'r');
if fid==-1
    fprintf(['\tFile NOT found for array positions at:\n\t' file_path '\configuration\config.txt' '!\n']);
else
    ARRAY_MICS = fscanf(fid,'%*s%*d%f%f',[2,64])';
    fclose(fid);
end

% Remove if necessary
if ~isempty(remove_mic)
    fprintf('\tRemoving defective mic data...\n');
    
    DATA(:, remove_mic) = [];
    ARRAY_MICS(remove_mic, :) = [];
    
else
    fprintf('\tNo defective mics removed!\n');
end

fprintf('\tRead in procedure completed.\n');
fprintf('\t------------------------------------------\n');

end