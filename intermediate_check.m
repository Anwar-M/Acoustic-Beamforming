%% Script for debugging CLEAN
close;
dynamic_range = 12; % dB
xmin = -.25; xmax = .25;
ymin = -.25; ymax = .25;
xmin = -1.5; xmax = 1.5;
ymin = -2; ymax = .5;

Pimp = P;
Pimp(real(Pimp)<0) = 0;

SPLint = 20*log10(sqrt(real(reshape(Pimp, N_X, N_Y).'))/2e-5);
maxval = ceil(max(real(SPLint(:))));
minval = maxval - dynamic_range;

figure('pos', [400 400 1500 500]);
subplot(1,3,1);SPLint(SPLint < minval - .5) = 0;
imagesc(X, Y, SPLint, [minval-.5 round(maxval)]);
title([num2str(norm(P)) ', ' num2str(norm(D))]);

set(0,'defaulttextinterpreter','latex');
hXLabel = xlabel('$x$ [m]');
hYLabel = ylabel('$y$ [m]');
set([hXLabel, hYLabel], ...
    'FontName', 'AvantGarde', 'FontSize', 14);

set(gca, ...
'YDir','normal', ...
'XTick', linspace(xmin, xmax, 5), ...
'YTick', linspace(ymin, ymax, 5), ...
'XColor'      , [.3 .3 .3], ...
'YColor'      , [.3 .3 .3], ...
'Fontsize'    , 16, ...
'FontName'   , 'Helvetica' );
set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', ':', 'xcolor', 'k', 'ycolor', 'k');

axis equal; axis([xmin xmax ymin ymax]);
cb = colorbar;
jetmod = jet;
jetmod(1,:) = 1;
colormap(jetmod);
cb.Limits = [minval maxval];

% Pomp = sum(conj(h).*(Cource*h), 1);
Pomp = sum(h.*(Cource*conj(h)), 1);
% Pomp(real(Pomp)<0) = 0;
SPLrem = 20*log10(sqrt(real(reshape(Pomp, N_X, N_Y).'))/2e-5);
maxval = ceil(max(real(SPLrem(:))));
minval = maxval - dynamic_range;

subplot(1,3,2);SPLrem(SPLrem < minval - .5) = 0;
imagesc(X, Y, SPLrem, [minval-.5 round(maxval)]);

title(num2str(norm(Q)));
set(0,'defaulttextinterpreter','latex');
hXLabel = xlabel('$x$ [m]');
hYLabel = ylabel('$y$ [m]');
set([hXLabel, hYLabel], ...
    'FontName', 'AvantGarde', 'FontSize', 14);

set(gca, ...
'YDir','normal', ...
'XTick', linspace(xmin, xmax, 5), ...
'YTick', linspace(ymin, ymax, 5), ...
'XColor'      , [.3 .3 .3], ...
'YColor'      , [.3 .3 .3], ...
'Fontsize'    , 16, ...
'FontName'   , 'Helvetica' );
set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', ':', 'xcolor', 'k', 'ycolor', 'k');

axis equal; axis([xmin xmax ymin ymax]);
cb = colorbar;
jetmod = jet;
jetmod(1,:) = 1;
colormap(jetmod);
cb.Limits = [minval maxval];

SPLflint = 20*log10(sqrt(real(reshape(Q, N_X, N_Y).'))/2e-5);
maxval = ceil(max(real(SPLflint(:))));
minval = maxval - dynamic_range;

subplot(1,3,3);SPLflint(SPLflint < minval - .5) = 0;
imagesc(X, Y, SPLflint, [minval-.5 round(maxval)]);

title(num2str(norm(Q)));
set(0,'defaulttextinterpreter','latex');
hXLabel = xlabel('$x$ [m]');
hYLabel = ylabel('$y$ [m]');
set([hXLabel, hYLabel], ...
    'FontName', 'AvantGarde', 'FontSize', 14);

set(gca, ...
'YDir','normal', ...
'XTick', linspace(xmin, xmax, 5), ...
'YTick', linspace(ymin, ymax, 5), ...
'XColor'      , [.3 .3 .3], ...
'YColor'      , [.3 .3 .3], ...
'Fontsize'    , 16, ...
'FontName'   , 'Helvetica' );
set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', ':', 'xcolor', 'k', 'ycolor', 'k');

axis equal; axis([xmin xmax ymin ymax]);
cb = colorbar;
jetmod = jet;
jetmod(1,:) = 1;
colormap(jetmod);
cb.Limits = [minval maxval];