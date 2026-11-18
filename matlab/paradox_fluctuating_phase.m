function paradox_fluctuating_phase()
% PARADOX_FLUCTUATING_PHASE  Section IV.E.
%   Jeltsema fluctuating-phase example extended to a balanced three-phase
%   four-wire load and identified with parallel-RL and parallel-RC windowed
%   local regressors.
%
%   Per phase:
%       v_x = sqrt(2) V cos(w0 t - theta_x)
%       i_x = sqrt(2) I cos(w0 t - theta_x - phi(t))
%       phi(t) = -gamma sin(w1 t),   w1 = w0/5.
%
%   Parameters used here follow Fig. 6 of Jeltsema et al.:
%       V = 1 V rms, I = 0.5 A rms, f0 = 50 Hz, gamma = pi/2.
%
%   Analytical single-phase formulas give the reference curves.  The
%   finite-window estimates solve the same local regressors:
%       parallel RL: i ~= G v + Gamma integral(v)
%       parallel RC: i ~= G v + C dv/dt.

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    p2_figstyle();
    set(groot, 'defaultAxesFontSize', 18);
    set(groot, 'defaultTextFontSize', 18);
    set(groot, 'defaultLegendFontSize', 14);
    set(groot, 'defaultLineLineWidth', 2.4);

    figdir = fullfile(here, '..', 'results', 'figures');
    resdir = fullfile(here, '..', 'results');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    if ~exist(resdir,'dir'); mkdir(resdir); end

    % Jeltsema Fig. 6 parameters.
    f0 = 50;
    om = 2*pi*f0;
    f1 = f0/5;
    om1 = 2*pi*f1;
    Vrms = 1;
    Irms = 0.5;
    Vm = sqrt(2)*Vrms;
    Im = sqrt(2)*Irms;
    gamma_amp = pi/2;

    fs = 200e3;
    Tend = 0.2;                         % two modulation periods
    t = (0:round(Tend*fs)-1).' / fs;    % endpoint excluded
    th = [0; -2*pi/3; +2*pi/3];

    va = Vm*cos(om*t - th(1));
    vb = Vm*cos(om*t - th(2));
    vc = Vm*cos(om*t - th(3));
    dv = [-Vm*om*sin(om*t - th(1)), ...
          -Vm*om*sin(om*t - th(2)), ...
          -Vm*om*sin(om*t - th(3))];
    qv = [Vm/om*sin(om*t - th(1)), ...
          Vm/om*sin(om*t - th(2)), ...
          Vm/om*sin(om*t - th(3))];

    phi = -gamma_amp*sin(om1*t);
    ia = Im*cos(om*t - th(1) - phi);
    ib = Im*cos(om*t - th(2) - phi);
    ic = Im*cos(om*t - th(3) - phi);

    vabc = [va, vb, vc];
    iabc = [ia, ib, ic];

    % Analytical reference from the instantaneous single-phase formulas.
    [G_RL, Gamma_ref] = analytic_parallel_rl(t, om, om1, Vrms, Irms, gamma_amp);
    [G_RC, C_ref] = analytic_parallel_rc(t, om, om1, Vrms, Irms, gamma_amp);

    % Finite-window estimates per phase. The plotted curves use phase a
    % because the analytical reference below is written for theta_a = 0.
    windowCyclesList = [0.10, 0.25, 0.5];
    hopCycles = 0.05;
    winRL = windowed_three_phase(vabc, iabc, qv, f0, fs, t, windowCyclesList, hopCycles);
    winRC = windowed_three_phase(vabc, iabc, dv, f0, fs, t, windowCyclesList, hopCycles);

    % Per-phase storage rate.  Its aggregate vanishes by balanced symmetry,
    % even though each phase exchanges stored energy.
    storage_current = Im*sin(phi);
    Wdot_a = va .* storage_current .* sin(om*t - th(1));
    Wdot_b = vb .* storage_current .* sin(om*t - th(2));
    Wdot_c = vc .* storage_current .* sin(om*t - th(3));
    Wdot_sum = Wdot_a + Wdot_b + Wdot_c;

    % Akagi p-q for aggregate context.
    K = sqrt(2/3);
    valpha = K*(va - 0.5*vb - 0.5*vc);
    vbeta  = K*(sqrt(3)/2*vb - sqrt(3)/2*vc);
    ialpha = K*(ia - 0.5*ib - 0.5*ic);
    ibeta  = K*(sqrt(3)/2*ib - sqrt(3)/2*ic);
    p_pq = valpha.*ialpha + vbeta.*ibeta;
    q_pq = valpha.*ibeta - vbeta.*ialpha;

    metrics = collect_metrics(G_RL, G_RC, Gamma_ref, C_ref, Wdot_a, Wdot_b, ...
        Wdot_c, Wdot_sum, p_pq, q_pq, winRL, winRC);

    fprintf('[Fluctuating phase] Jeltsema three-phase extension: V=%.3g Vrms, I=%.3g Arms, f0=%.1f Hz, f1=%.1f Hz, gamma=pi/2\n', ...
        Vrms, Irms, f0, f1);
    fprintf('[Fluctuating phase] analytic G range [%.4g, %.4g] S; Gamma [%.2f, %.2f] 1/H; C [%.4g, %.4g] F\n', ...
        metrics.G_range, metrics.Gamma_range, metrics.C_range);
    fprintf('[Fluctuating phase] rms(sum Wdot_x)=%.3e W; phase rms Wdot=[%.3f %.3f %.3f] W\n', ...
        metrics.Wdot_sum_rms, metrics.Wdot_phase_rms);
    fprintf('[Fluctuating phase] mean(p_pq)=%.4f W, mean(q_pq)=%.3e var, rms(q_pq)=%.4f var\n', ...
        mean(p_pq), mean(q_pq), rms_local(q_pq));
    for kidx = 1:numel(windowCyclesList)
        fprintf('[Fluctuating phase] %.2f-cycle windows: RL residual %.3e, RC residual %.3e\n', ...
            windowCyclesList(kidx), median(winRL(kidx).residual(:), 'omitnan'), ...
            median(winRC(kidx).residual(:), 'omitnan'));
    end

    fh = figure('Position',[80 80 720 620],'Color','w');
    layout = tiledlayout(fh, 4, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % Row 1: terminal context — phase a waveforms (left) and modulation phi (right)
    nexttile(layout);
    yyaxis left;  plot(t*1e3, va); ylabel('$v_a$ [V]'); set(gca,'YColor','k');
    yyaxis right; plot(t*1e3, ia); ylabel('$i_a$ [A]'); set(gca,'YColor','k');
    title(sprintf('Phase $a$ waveforms ($V{=}%.1f$~V rms, $I{=}%.1f$~A rms)', Vrms, Irms));
    grid on;

    nexttile(layout);
    plot(t*1e3, phi*180/pi);
    ylabel('$\varphi(t)$ [deg]');
    title('Modulation $\varphi(t){=}-\frac{\pi}{2}\sin\omega_1 t$, $\omega_1{=}\omega_0/5$');
    grid on; ylim_pad(gca, 0.08);

    % Row 2: identified G(t) under each topology — must coincide
    nexttile(layout);
    plot(t*1e3, G_RL, 'k-'); hold on;
    plot_window_series(winRL, 1, 1, 'G', 1);
    ylabel('$G_a(t)$ [S]');
    title('Parallel $RL$: analytical $G_a(t)$ and windowed estimate');
    legend({'analytic', 'windowed'}, 'Location', 'northeast');
    grid on; ylim_pad(gca, 0.08);

    nexttile(layout);
    plot(t*1e3, G_RC, 'k-'); hold on;
    plot_window_series(winRC, 1, 1, 'G', 1);
    ylabel('$G_a(t)$ [S]');
    title('Parallel $RC$: same $G_a(t)$ recovered');
    legend({'analytic', 'windowed'}, 'Location', 'northeast');
    grid on; ylim_pad(gca, 0.08);

    % Row 3: storage coordinate under each topology — mirror images
    nexttile(layout);
    plot(t*1e3, Gamma_ref, 'k-'); hold on;
    plot_window_series(winRL, 1, 1, 'storage', 1);
    ylabel('$\Gamma_a(t)=1/L_a$ [1/H]');
    title('Parallel $RL$ storage coordinate');
    legend({'analytic', 'windowed'}, 'Location', 'northeast');
    grid on; ylim_pad(gca, 0.12);

    nexttile(layout);
    plot(t*1e3, C_ref*1e3, 'k-'); hold on;
    plot_window_series(winRC, 1, 1, 'storage', 1e3);
    ylabel('$C_a(t)$ [mF]');
    title('Parallel $RC$ storage coordinate ($C{=}-\Gamma/\omega_0^2$)');
    legend({'analytic', 'windowed'}, 'Location', 'northeast');
    grid on; ylim_pad(gca, 0.12);

    % Row 4: per-phase storage rates (spanned)
    nexttile(layout, [1 2]);
    plot(t*1e3, Wdot_a); hold on;
    plot(t*1e3, Wdot_b);
    plot(t*1e3, Wdot_c);
    xlabel('time [ms]'); ylabel('$W_x^{\prime}$ [W]');
    title(sprintf('Per-phase storage rates ($\\sum_x W_x^{\\prime}{\\equiv}0$ in balanced 3-phase, rms ${=}%.1e$ W)', ...
        metrics.Wdot_sum_rms));
    legend({'$W_a^{\prime}$', '$W_b^{\prime}$', '$W_c^{\prime}$'}, ...
        'Location', 'east');
    grid on; ylim_pad(gca, 0.08);

    out_pdf = fullfile(figdir,'fig_p2_fluctuating_phase.pdf');
    exportgraphics(fh, out_pdf,'ContentType','vector');
    close(fh);

    save(fullfile(resdir,'paradox_fluctuating_phase.mat'), ...
        't','va','vb','vc','ia','ib','ic','phi', ...
        'G_RL','G_RC','Gamma_ref','C_ref','winRL','winRC', ...
        'Wdot_a','Wdot_b','Wdot_c','Wdot_sum', ...
        'p_pq','q_pq','metrics','om','om1','gamma_amp','Vrms','Irms');

    fprintf('[Fluctuating phase] figure: %s\n', out_pdf);
end

function [G, Gamma] = analytic_parallel_rl(t, om, om1, Vrms, Irms, gamma_amp)
    Vm = sqrt(2)*Vrms;
    Im = sqrt(2)*Irms;
    phi = -gamma_amp*sin(om1*t);
    psi = om*t - phi;
    dphi = -gamma_amp*om1*cos(om1*t);
    dpsi = om - dphi;

    v = Vm*cos(om*t);
    qv = Vm/om*sin(om*t);
    dv = -Vm*om*sin(om*t);
    i = Im*cos(psi);
    di = -Im*sin(psi).*dpsi;

    den = v.^2 - qv.*dv;
    G = (i.*v - qv.*di)./den;
    Gamma = (v.*di - i.*dv)./den;
end

function [G, C] = analytic_parallel_rc(t, om, om1, Vrms, Irms, gamma_amp)
    Vm = sqrt(2)*Vrms;
    Im = sqrt(2)*Irms;
    phi = -gamma_amp*sin(om1*t);
    psi = om*t - phi;
    dphi = -gamma_amp*om1*cos(om1*t);
    dpsi = om - dphi;

    v = Vm*cos(om*t);
    dv = -Vm*om*sin(om*t);
    ddv = -Vm*om^2*cos(om*t);
    i = Im*cos(psi);
    di = -Im*sin(psi).*dpsi;

    den = v.*ddv - dv.^2;
    G = (i.*ddv - dv.*di)./den;
    C = (v.*di - i.*dv)./den;
end

function win = windowed_three_phase(vabc, iabc, x2abc, f0, fs, t, windowCyclesList, hopCycles)
    win = repmat(struct('windowCycles', [], 'time', [], 'G', [], ...
        'storage', [], 'residual', []), numel(windowCyclesList), 1);
    for k = 1:numel(windowCyclesList)
        windowSamples = max(12, round(windowCyclesList(k)*fs/f0));
        hopSamples = max(1, round(hopCycles*fs/f0));
        starts = 1:hopSamples:(numel(t)-windowSamples+1);
        G = NaN(numel(starts), 3);
        storage = NaN(numel(starts), 3);
        residual = NaN(numel(starts), 3);
        center = NaN(numel(starts), 1);
        for w = 1:numel(starts)
            idx = starts(w):(starts(w)+windowSamples-1);
            center(w) = mean(t(idx));
            for ph = 1:3
                X = [vabc(idx,ph), x2abc(idx,ph)];
                y = iabc(idx,ph);
                theta = normalized_solve(X, y);
                G(w,ph) = theta(1);
                storage(w,ph) = theta(2);
                residual(w,ph) = norm(y - X*theta)/max(norm(y), eps);
            end
        end
        win(k).windowCycles = windowCyclesList(k);
        win(k).time = center;
        win(k).G = G;
        win(k).storage = storage;
        win(k).residual = residual;
    end
end

function theta = normalized_solve(X, y)
    colScale = sqrt(mean(X.^2, 1)).';
    yScale = sqrt(mean(y.^2));
    colScale(~isfinite(colScale) | colScale <= eps) = 1;
    yScale = max(yScale, eps);
    thetaScaled = (X./colScale.') \ (y./yScale);
    theta = thetaScaled .* yScale ./ colScale;
end

function metrics = collect_metrics(G_RL, G_RC, Gamma_ref, C_ref, Wdot_a, Wdot_b, ...
    Wdot_c, Wdot_sum, p_pq, q_pq, winRL, winRC)
    metrics = struct();
    metrics.G_range = [min([G_RL; G_RC]), max([G_RL; G_RC])];
    metrics.G_cross_rms = rms_local(G_RL - G_RC);
    metrics.Gamma_range = [min(Gamma_ref), max(Gamma_ref)];
    metrics.C_range = [min(C_ref), max(C_ref)];
    metrics.Wdot_phase_rms = [rms_local(Wdot_a), rms_local(Wdot_b), rms_local(Wdot_c)];
    metrics.Wdot_sum_rms = rms_local(Wdot_sum);
    metrics.p_mean = mean(p_pq);
    metrics.q_mean = mean(q_pq);
    metrics.q_rms = rms_local(q_pq);
    metrics.windowCycles = [winRL.windowCycles];
    metrics.RL_residual = arrayfun(@(s)median(s.residual(:),'omitnan'), winRL);
    metrics.RC_residual = arrayfun(@(s)median(s.residual(:),'omitnan'), winRC);
end

function plot_window_series(win, idx, phaseIdx, field, scale)
    values = win(idx).(field)(:,phaseIdx) * scale;
    plot(win(idx).time*1e3, values, 'o', 'MarkerSize', 4, 'LineWidth', 1.1);
end

function y = rms_local(x)
    y = sqrt(mean(x.^2, 'omitnan'));
end

function ylim_pad(ax, frac)
    yl = ylim(ax);
    span = yl(2) - yl(1);
    if span <= 0
        span = max(abs(yl), [], 'all');
    end
    pad = max(frac*span, eps);
    ylim(ax, [yl(1)-pad, yl(2)+pad]);
end
