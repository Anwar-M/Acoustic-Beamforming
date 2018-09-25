function plot_settings_font(axisHandle, ...
                       xlab, ylab, titlab, ...
                       xlim, ylim, xtick, ytick, fontsize, ...
                       xgrid, ygrid, ax_equal, clr_bar, ctitlab,...
                       save_imgs, filepath)
% Customize your plots 'easily'.
%
% Anwar Malgoezar, Jan. 2017
% Group ANCE

% Needs CMU Serif font for consistency

% x, y and title labels. Paper-ready
hXLabel = xlabel(axisHandle, xlab, 'Interpreter', 'LaTex');
hYLabel = ylabel(axisHandle, ylab, 'Interpreter', 'LaTex');
hTitle = title(axisHandle, titlab, 'Interpreter', 'LaTex');

% Equal axis, or not
if ax_equal axis(axisHandle, 'equal'); end;

set([hXLabel, hYLabel, hTitle], 'FontSize', fontsize);

set(axisHandle, ...
      'FontName'    , 'CMU Serif', ...
      'YDir'        , 'normal'   , ...
      'Box'         , 'on'       , ...
      'TickDir'     , 'in'       , ...
      'TickLength'  , [.01 .01]  , ...
      'XGrid'       , xgrid      , ...
      'YGrid'       , ygrid      , ...
      'XColor'      , [.3 .3 .3] , ...
      'YColor'      , [.3 .3 .3] , ...
      'Fontsize'    , fontsize         , ...
      'LineWidth'   , 1          , ...
	  'TickLabelInterpreter', 'latex');

if (length(xlim)==2)&&(length(ylim)==2)
    set(axisHandle, ...
      'XLim'        , [xlim(1) xlim(2)], ...
      'YLim'        , [ylim(1) ylim(2)], ...
      'XTick'       , xtick, ...
      'YTick'       , ytick);
end

if clr_bar(1)
    %cb = colorbar;
    cb = colorbar('YTick', linspace(clr_bar(2),clr_bar(3),7));
    caxis([clr_bar(2) clr_bar(3)])
    title(cb, ctitlab, 'Fontsize', fontsize, 'FontName', 'CMU Serif', 'Interpreter', 'Latex');
%     ylabel(cb, ctitlab, 'Fontsize', fontsize, 'FontName', 'CMU Serif', 'Interpreter', 'Latex');
    set(cb, 'TickLabelInterpreter', 'latex');
    
end

% scalefactor = .85;
% g = get(gca,'Position');
% g(1:2) = g(1:2) + (1-scalefactor)/2*g(3:4);
% g(3:4) = scalefactor*g(3:4);
% set(gca,'Position',g);
    
if save_imgs
    set(gcf, 'PaperUnits','Points')
%     pos = get(gcf, 'Position');
%     set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)]);
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf,'renderer','opengl');
    print(filepath, '-dpng', '-r300');
    print(filepath, '-depsc', '-r300');
    hgsave(filepath);
end

end