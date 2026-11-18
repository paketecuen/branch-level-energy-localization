function paradox_balanced_RL_compare()
% PARADOX_BALANCED_RL_COMPARE  §IV.A enhancement.
%   Phase a of a balanced 3-phase RL load.
%   Branch-level decomposition (energy-balance) vs IEEE Std. 1459
%   instantaneous active and reactive per phase. Same amplitudes, with a
%   phase shift of 2*theta in the active terms and of pi - 2*theta in the
%   reactive terms, where theta = atan(omega L / R).

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    p2_figstyle();
    % Local style overrides (Fig 3 only): heavier lines + larger labels
    set(groot, 'defaultAxesFontSize', 20);
    set(groot, 'defaultTextFontSize', 20);
    set(groot, 'defaultLegendFontSize', 18);
    set(groot, 'defaultLineLineWidth', 2.6);
    figdir = fullfile(here, '..', 'results', 'figures');
    resdir = fullfile(here, '..', 'results');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    if ~exist(resdir,'dir'); mkdir(resdir); end

    f0 = 50;
    om = 2*pi*f0;
    fs = 50e3;
    Tend = 3/f0;
    t  = (0:1/fs:Tend).';

    V = 120;
    R = 1.1;
    L = 1e-3;
    Z = sqrt(R^2 + (om*L)^2);
    I = V/Z;
    th_z = atan(om*L/R);

    P1 = V*I*cos(th_z);
    Q1 = V*I*sin(th_z);

    % Phase a
    va = sqrt(2)*V*sin(om*t);
    ia = sqrt(2)*I*sin(om*t - th_z);
    iad = sqrt(2)*I*om*cos(om*t - th_z);

    % Branch-level decomposition (energy-balance per branch)
    p_branch_act = R*ia.^2;
    p_branch_rea = L*ia.*iad;

    % IEEE Std. 1459 instantaneous per-phase formulas (Section 3.1.1.1)
    P = P1; Q = Q1;
    p_1459_act = P*(1 - cos(2*om*t));
    p_1459_rea = -Q*sin(2*om*t);

    % Sanity: both pairs sum to v_a*i_a
    p_term = va.*ia;
    err_branch = norm(p_term - (p_branch_act + p_branch_rea))/norm(p_term);
    err_1459   = norm(p_term - (p_1459_act + p_1459_rea))/norm(p_term);

    fprintf('[Balanced RL compare] theta = %.2f deg, P1 = %.2f W, Q1 = %.2f var\n', ...
        th_z*180/pi, P1, Q1);
    fprintf('[Balanced RL compare] sum residual branch: %.3e ; 1459: %.3e\n', err_branch, err_1459);

    fh = figure('Position',[100 100 680 490],'Color','w');

    subplot(2,1,1);
    plot(t*1e3, p_branch_act/1e3,'-'); hold on;
    plot(t*1e3, p_1459_act/1e3,'--');
    legend({'$p_{R,a}=R\,i_a^2$', ...
            '$P[1-\cos 2\omega t]$'}, ...
        'Location','northeast');
    ylabel('active [kW]');
    title(sprintf('Phase $a$ active power: same amplitude, phase shift $2\\theta\\approx %.0f^\\circ$', ...
        2*th_z*180/pi));
    grid on;

    subplot(2,1,2);
    plot(t*1e3, p_branch_rea/1e3,'-'); hold on;
    plot(t*1e3, p_1459_rea/1e3,'--');
    legend({'$W_{L,a}^{\prime}=L\,i_a i_a^{\prime}$', ...
            '$-Q\sin 2\omega t$'}, ...
        'Location','northeast');
    xlabel('time [ms]'); ylabel('reactive [kvar]');
    title(sprintf('Phase $a$ reactive component: same amplitude, phase shift $\\pi-2\\theta\\approx %.0f^\\circ$', ...
        180 - 2*th_z*180/pi));
    grid on;

    out_pdf = fullfile(figdir,'fig_p2_balanced_RL_compare.pdf');
    exportgraphics(fh, out_pdf,'ContentType','vector');
    close(fh);

    save(fullfile(resdir,'paradox_balanced_RL_compare.mat'), ...
        't','va','ia','p_branch_act','p_branch_rea','p_1459_act','p_1459_rea', ...
        'V','R','L','I','th_z','P1','Q1');

    fprintf('[Balanced RL compare] figure: %s\n', out_pdf);
end
