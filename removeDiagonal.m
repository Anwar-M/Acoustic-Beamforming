function CSM = removeDiagonal(CSM)
% Removes the diagonal from the CSM. For info see:
% "Detection of aeroacoustic sound sources on aircraft and wind turbines" 
% By S. Oerlemans
%
% Assumes a CSM as a 2D or 3D matrix of n_mic x n_mic x n_f
% n_mic: amount of microphones
% n_f: amount of frequencies, can be 1, i.e. a 2D matrix
%
% Remember to remove (i.e. = 0, or better < 1e-20 to prevent -Inf dB values) 
% negative pressures from beamform map! The beamform map are the source
% powers obtained for you scan plane immediately after beamforming and
% before calculating SPL values.
%
% Anwar Malgoezar, June 2017
% Group ANCE

n_mic = size(CSM,1);
n_f = size(CSM,3);
CSMvec = CSM(:);
for ff = 1:n_f
    CSMvec((n_mic+1)*(0:(n_mic-1))+1+n_mic*n_mic*(ff-1)) = 0;
end
CSM = reshape(CSMvec,n_mic,n_mic,n_f);
clear CSMvec;