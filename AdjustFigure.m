figfile = 'O:\V-Tunnel 13-12 Thrust Experiment\results 15-12\spectra_propeller_static_85.fig';
anotherfig = 'O:\V-Tunnel 21-07 Plate Duct Shielding\Results\spectra_monopole.fig';
save_data = 1;

close all;
addpath('O:\MATLAB Signal Processing Files');
open(anotherfig);


h = findobj(gca,'Type','line');
X = h(3).XData;
Y = h(3).YData;
% h(1).LineWidth = 1;
% h(2).LineWidth = 1;
% h(3).LineWidth = 1;
% h(3).Color = [0, 0.4470, 0.7410];
% h(2).Color = [0.8500, 0.3250, 0.0980];
% h(1).Color = [0.4660, 0.6740, 0.1880];
% h(6).Color = [0, 0.4470, 0.7410];
% h(5).Color = [0.8500, 0.3250, 0.0980];
% h(4).Color = [0.4660, 0.6740, 0.1880];
% h(4).LineWidth = 1;
% h(5).LineWidth = 1;
% h(6).LineWidth = 1;
close;
open(figfile);
hold on;
plot(X, Y, 'LineWidth', 1, 'Color', [1 1 1]*0.7);
hold off;

% s = findobj('type','legend');
% delete(s);
% 
% hl = legend('Duct', 'No Duct', 'Background', 'Bla');
% set(hl, 'Interpreter', 'LaTex');

plot_settings(gca, '$f$ [kHz]', 'PSD [dB/Hz]', [], [0 5], [-2 85], ...
    0:1e0:5e0, 0:10:80, 'on', 'on', 0, 0, [], save_data, ...
    [figfile(1:end-4) 'mod']);

%%
diro = 'O:\PhD Thesis\RESULTS\CH5\LOUDSPEAKER';
listt = dir([diro '\BF*.fig']);
for I = 1:length(listt)
    open([diro '\' listt(I).name]);
    set(gcf, 'PaperUnits','Points')
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf,'renderer','opengl');
    print([diro '\' listt(I).name(1:end-4)], '-depsc', '-r300');
    close;
end

%%
diro = 'O:\PhD Thesis\RESULTS\CH5\SPARKER';
listt = dir([diro '\BF*noflow*.fig']);
for I = 1:length(listt)
    open([diro '\' listt(I).name]);
    h = findobj(gca);
    h(6).LevelList = [h(6).LevelList(1):h(6).LevelList(end)];
    set(gcf, 'PaperUnits','Points')
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf,'renderer','opengl');
    print([diro '\' listt(I).name(1:end-4)], '-depsc', '-r300');
    close;
end
