t = 16;

%%
tic;
segsize = t*50e3;
file = 'acoustic_data';
fid = fopen(file, 'r');
STATUS = fseek(fid, 384000000, 0);
curr_data = fread(fid, [64 segsize], 'int16', 0, 'l');
fclose(fid);
toc

%%

fid = fopen('acoustic_data_trimmed', 'wb');
fwrite(fid, curr_data, 'int16');
fclose(fid);