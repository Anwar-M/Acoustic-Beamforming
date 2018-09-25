function [CSM, freqs] = developCSM(x, freq_l, freq_u, Fs, t_block, overlap, t_start, t_end)
% Develop CSM from microphone array measurements. Pressure time signals.
% Creates time blocks and averages to reduce noise for CSM given your
% frequency range with possible overlap to generate more blocks. See
% Welch's method for more info.
%
% Anwar Malgoezar, June 2015
% Group ANCE

% Time series as columns!
fprintf('\t------------------------------------------\n');
fprintf('\tStart building CSM...\n');

N_total_samples = size(x, 1);
N_signals = size(x, 2);
t_signal = N_total_samples/Fs;

fprintf('\tSignal time: %.2f s, Sample frequency: %d Hz,\n\t%d samples and %d signals...\n', t_signal, Fs, N_total_samples, N_signals);

if nargin < 8
    t_start = 0;
    t_end = t_signal;
end
if (t_start < 0) || (t_end > t_signal) || (t_block>(t_end-t_start))
    error('Time-boundaries out of bounds!');
end

fprintf('\tUsing signal between %.2f and %.2f s...\n', t_start, t_end);
fprintf('\tBlock-size %.2f s with %.2f%% overlap.\n', t_block, overlap*100);

start_sample = floor(t_start*Fs) + 1;
block_samples = ceil(t_block*Fs);
end_sample = ceil(t_end*Fs);

% Determine relative OFFSET sample value depending on OVERLAP criterion.
% I.e. if we have 5000 sample block with 50% overlap, the next START sample
% for the NEXT block would be at sample 5000/2 + 1 = 2501, so at half of 
% the current block.
offset_sample = block_samples * (1 - overlap) + 1 - 1;

% Calculate amount of blocks able to fit between start and end time with
% overlap
N_blocks = floor( (end_sample - block_samples - start_sample + 1) / ...
                  offset_sample ) + 1;
              
x_fr = Fs / block_samples * (0:floor(block_samples/2)-1);
freq_sels = find((x_fr>=freq_l).*(x_fr<=freq_u));
N_freqs = length(freq_sels);

fprintf('\tFrequency range: %.2f to %.2f Hz, %d frequency points...\n', freq_l, freq_u, N_freqs);
CSM = zeros(N_signals, N_signals, N_freqs);

reverseStr = '';
for B = 1:N_blocks
    msg = sprintf('\tEvaluating CSM at block %d/%d...\n', B, N_blocks);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    N_start = start_sample + (B-1) * offset_sample;
    N_end = N_start + block_samples - 1;
    x_fft = 2*fft(x(N_start:N_end,:))/block_samples;
%     x_fft = 2*fft(x(N_start:N_end,:).*(hann(numel(N_start:N_end))*ones(1,size(x,2))))/block_samples;
    for F = 1:N_freqs
        CSM(:,:,F) = CSM(:,:,F) + 0.5*x_fft(freq_sels(F),:)'*x_fft(freq_sels(F),:);
    end
    
end

CSM = CSM/N_blocks;
freqs = x_fr(freq_sels);
fprintf('\tFinished building CSM!\n');
fprintf('\t------------------------------------------\n');

end