function p2_figstyle()
% P2_FIGSTYLE  Apply consistent LaTeX + large font style to current MATLAB session.
%   Call once at the start of any plotting script.

    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

    set(groot, 'defaultAxesFontName', 'Times');
    set(groot, 'defaultTextFontName', 'Times');

    set(groot, 'defaultAxesFontSize', 20);
    set(groot, 'defaultTextFontSize', 20);
    set(groot, 'defaultLegendFontSize', 18);
    set(groot, 'defaultColorbarFontSize', 18);

    set(groot, 'defaultLineLineWidth', 1.6);
    set(groot, 'defaultAxesLineWidth', 1.0);
    set(groot, 'defaultAxesBox', 'on');
    set(groot, 'defaultAxesGridAlpha', 0.25);
end
