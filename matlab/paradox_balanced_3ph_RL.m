function paradox_balanced_3ph_RL()
% PARADOX_BALANCED_3PH_RL  Section IV.A: Balanced three-phase RL load.
%   Aggregated terminal power and aggregated magnetic energy are constant
%   in steady state, while per-phase storage rates exchange cyclically.

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    p2_figstyle();
    % Local style overrides (Fig 2 only): heavier lines + larger labels
    set(groot, 'defaultAxesFontSize', 20);
    set(groot, 'defaultTextFontSize', 20);
    set(groot, 'defaultLegendFontSize', 18);
    set(groot, 'defaultLineLineWidth', 2.6);
    figdir = fullfile(here, '..', 'results', 'figures');
    resdir = fullfile(here, '..', 'results');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    if ~exist(resdir,'dir'); mkdir(resdir); end

    % System parameters
    f0 = 50;
    om = 2*pi*f0;
    Vm = 230*sqrt(2);
    R  = 5.0;
    L  = 0.020;
    Z  = R + 1j*om*L;
    Im = Vm/abs(Z);
    phi = angle(Z);

    fs = 20e3;
    t  = (0:1/fs:3/f0).';

    th = [0; -2*pi/3; +2*pi/3];

    % Phase voltages and currents
    va = Vm*cos(om*t - th(1));
    vb = Vm*cos(om*t - th(2));
    vc = Vm*cos(om*t - th(3));
    ia = Im*cos(om*t - th(1) - phi);
    ib = Im*cos(om*t - th(2) - phi);
    ic = Im*cos(om*t - th(3) - phi);

    % Per-phase analytic derivatives
    iad = -Im*om*sin(om*t - th(1) - phi);
    ibd = -Im*om*sin(om*t - th(2) - phi);
    icd = -Im*om*sin(om*t - th(3) - phi);

    % Branch energy terms
    pR_a = R*ia.^2;  pR_b = R*ib.^2;  pR_c = R*ic.^2;
    pR_tot = pR_a + pR_b + pR_c;

    Wdot_a = L*ia.*iad;
    Wdot_b = L*ib.*ibd;
    Wdot_c = L*ic.*icd;
    Wdot_tot = Wdot_a + Wdot_b + Wdot_c;

    WL_a = 0.5*L*ia.^2;
    WL_b = 0.5*L*ib.^2;
    WL_c = 0.5*L*ic.^2;
    WL_tot = WL_a + WL_b + WL_c;

    p_term = va.*ia + vb.*ib + vc.*ic;

    fprintf('[Balanced 3ph RL] max ripple of total power: %.3e W\n', max(p_term)-min(p_term));
    fprintf('[Balanced 3ph RL] max ripple of total stored energy: %.3e J\n', max(WL_tot)-min(WL_tot));
    fprintf('[Balanced 3ph RL] max |sum dW/dt|: %.3e W\n', max(abs(Wdot_tot)));

    fh = figure('Position',[100 100 680 640],'Color','w');

    subplot(4,1,1);
    plot(t*1e3, va); hold on; plot(t*1e3, vb); plot(t*1e3, vc);
    legend({'$v_a$','$v_b$','$v_c$'}, 'Location','northeast');
    ylabel('voltage [V]'); title('Balanced phase voltages'); grid on;

    subplot(4,1,2);
    plot(t*1e3, p_term/1e3,'-'); hold on;
    plot(t*1e3, pR_tot/1e3,'--');
    legend({'$p_{3\phi}$','$\sum p_R$'}, 'Location','east');
    ylabel('[kW]'); ylim([0 15]);
    title('Aggregated terminal and Joule powers (constant)'); grid on;

    subplot(4,1,3);
    plot(t*1e3, WL_a*1e3); hold on;
    plot(t*1e3, WL_b*1e3); plot(t*1e3, WL_c*1e3);
    plot(t*1e3, WL_tot*1e3,'k--');
    legend({'$W_{L,a}$','$W_{L,b}$','$W_{L,c}$','$\sum W_L$'}, 'Location','east');
    ylabel('[mJ]'); title('Per-phase magnetic energies (oscillate; sum constant)'); grid on;

    subplot(4,1,4);
    plot(t*1e3, Wdot_a/1e3); hold on;
    plot(t*1e3, Wdot_b/1e3); plot(t*1e3, Wdot_c/1e3);
    plot(t*1e3, Wdot_tot/1e3,'k--');
    legend({'$W_{L,a}^{\prime}$','$W_{L,b}^{\prime}$','$W_{L,c}^{\prime}$','$\sum W_L^{\prime}$'}, 'Location','east');
    ylabel('[kW]'); xlabel('time [ms]');
    title('Per-phase storage rates exchange cyclically; aggregated $\sum W_L^{\prime}=0$');
    grid on;

    out_pdf = fullfile(figdir,'fig_p2_balanced_3ph_rl_energy.pdf');
    exportgraphics(fh, out_pdf,'ContentType','vector');
    close(fh);

    save(fullfile(resdir,'paradox_balanced_3ph_RL.mat'), ...
        't','va','vb','vc','ia','ib','ic','pR_tot','Wdot_a','Wdot_b','Wdot_c', ...
        'Wdot_tot','WL_a','WL_b','WL_c','WL_tot','p_term','R','L','om');

    fprintf('[Balanced 3ph RL] figure: %s\n', out_pdf);
end
