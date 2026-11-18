function fourwire_nonlinear()
% FOURWIRE_NONLINEAR  Section IV-F. Unbalanced four-wire nonlinear load.
%
%   Phase a: hysteretic inductor in parallel with a capacitor C_a.
%   Phase b: linear parallel R_b, L_b, C_b branch.
%   Phase c: series R_c-L_c branch with a freewheeling diode, connected
%            through a transistor that conducts for 7 ms in each period.
%   Supply:  v_x = 120*sqrt(2)*sin(w t - theta_x) + 0.12*sqrt(2)*sin(3 w t),
%            theta = (0, 2*pi/3, -2*pi/3), f = 50 Hz, neutral connected.
%
%   The waveforms are computed from this circuit specification. The script
%   evaluates the branch-level decomposition of each phase, the neutral
%   current, and the series-parallel map of Corollary 2 for phase b. The
%   result is compared with the curves of the two four-wire figures of the
%   manuscript (data/fourwire_reference), which were obtained in 2021 from
%   the parameters identified on the same circuit.
%
%   The series elastance S_s of Corollary 2 depends on the integration
%   constant of the primitive of the current. The zero-mean primitive is used
%   here, whereas the corresponding panel of the manuscript was computed with
%   a primitive that includes a large constant. R_s, L_s, and the three power
%   terms do not depend on that constant.

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    p2_figstyle();
    set(groot, 'defaultAxesFontSize', 14);
    set(groot, 'defaultTextFontSize', 14);
    set(groot, 'defaultLegendFontSize', 12);
    set(groot, 'defaultLineLineWidth', 1.4);
    root   = fullfile(here, '..');
    figdir = fullfile(root, 'results', 'figures');
    resdir = fullfile(root, 'results');
    refdir = fullfile(root, 'data', 'fourwire_reference');
    if ~exist(figdir,'dir'); mkdir(figdir); end

    % ---------------------------------------------------------- parameters
    par.f0 = 50;  par.w = 2*pi*par.f0;  par.T = 1/par.f0;
    par.Vm = 120*sqrt(2);  par.V3 = 0.12*sqrt(2);
    par.theta = [0, 2*pi/3, -2*pi/3];        % v_x = Vm*sin(w t - theta_x) + V3*sin(3 w t)
    % Phase a: capacitor and piecewise-linear hysteresis loop of the inductor
    par.Ca = 500e-6;                         % [F]
    par.Lsteep = 0.01;  par.Lflat = 0.001;   % slopes d(phi)/d(i) of the loop [H]
    par.phi1 = 0.15;    par.phi2 = 0.40;     % flux linkages of the knees [Wb]
    % Phase b: linear parallel RLC
    par.Rb = 2.0;  par.Lb = 0.01;  par.Cb = 2500e-6;
    % Phase c: switched series RL
    par.Rc = 0.75;  par.Lc = 0.4e-3;
    par.ton = 15.33e-3;  par.toff = 22.33e-3;   % conduction interval in each period [s]

    % ---------------------------------------------- waveforms on one period
    N = 10000;
    t = (0:N-1).' * par.T/N;
    A = phase_a(t, par);
    B = phase_b(t, par);
    C = phase_c(t, par);
    i_n = -(A.i + B.i + C.i);                % neutral current

    relres = @(s) norm(s.p - (s.pR + s.qL + s.qC)) / norm(s.p);
    fprintf('[4-wire] balance residual |p - (p_R + q_L + q_C)|/|p|: a = %.2e, b = %.2e, c = %.2e\n', ...
        relres(A), relres(B), relres(C));
    fprintf('[4-wire] average power P:   a = %8.1f W, b = %8.1f W, c = %8.1f W\n', mean(A.p), mean(B.p), mean(C.p));
    fprintf('[4-wire] average of p_R:    a = %8.1f W, b = %8.1f W, c = %8.1f W\n', mean(A.pR), mean(B.pR), mean(C.pR));
    fprintf('[4-wire] average of q_L:    a = %8.1f W, b = %8.1f W, c = %8.1f W\n', mean(A.qL), mean(B.qL), mean(C.qL));
    fprintf('[4-wire] average of q_C:    a = %8.1f W, b = %8.1f W, c = %8.1f W\n', mean(A.qC), mean(B.qC), mean(C.qC));
    Ehyst = sum(A.qL) * par.T/N;             % energy absorbed by the hysteretic inductor per cycle
    fprintf('[4-wire] phase a: hysteresis loop area = %.3f J per cycle, i.e. %.1f W at %g Hz (no resistor)\n', ...
        Ehyst, Ehyst*par.f0, par.f0);
    ioff = max(C.i(C.on & circshift(~C.on, -1)));   % current at turn-off
    fprintf('[4-wire] phase c: conduction %.2f-%.2f ms, peak current %.2f A, current at turn-off %.2f A\n', ...
        par.ton*1e3, par.toff*1e3, max(C.i), ioff);
    fprintf('[4-wire] phase c: energy left in L_c at turn-off = %.3f J per cycle (%.1f W), dissipated in R_c through the diode\n', ...
        0.5*par.Lc*ioff^2, 0.5*par.Lc*ioff^2*par.f0);
    fprintf('[4-wire] rms currents: a = %.1f A, b = %.1f A, c = %.1f A, neutral = %.1f A\n', ...
        rms_local(A.i), rms_local(B.i), rms_local(C.i), rms_local(i_n));

    % ------------------------- series-parallel map of Corollary 2, phase b
    D = duality_phase_b(B, par);
    fprintf('[4-wire] Corollary 2 on phase b: series equation residual max|R_s i + L_s i'' + S_s q - v|/max|v| = %.2e\n', D.kvl);
    fprintf('[4-wire] Corollary 2 on phase b: max |parallel - series| = %.2e W (p_R), %.2e W (q_L), %.2e W (q_C)\n', ...
        D.err_pR, D.err_qL, D.err_qC);
    fprintf('[4-wire] Corollary 2 on phase b: integral of L_s i i'' vs (1/2) Gamma phi^2: max difference %.2e J (peak energy %.2f J)\n', ...
        D.err_WL, D.WLmax);
    fprintf('[4-wire] Corollary 2 on phase b: R_s >= 0 everywhere: %d; L_s changes sign: %d; S_s changes sign: %d\n', ...
        D.Rs_nonneg, D.Ls_sign, D.Ss_sign);

    % ---------------------------- comparison with the curves of the figures
    cmp = compare_with_reference(refdir, par);

    % ---------------------------------------------------------------- figures
    plot_powers(par, fullfile(figdir, 'fig_fourwire_powers.pdf'));
    plot_duality(par, fullfile(figdir, 'fig_fourwire_duality.pdf'));

    save(fullfile(resdir, 'fourwire_nonlinear.mat'), 't', 'A', 'B', 'C', 'i_n', 'D', 'par', 'cmp');
    fprintf('[4-wire] figures: fig_fourwire_powers.pdf, fig_fourwire_duality.pdf\n');
end

% ======================================================================
function [v, dv, phi, ddv, phi2] = supply(t, th, par)
% Phase voltage, its derivatives and its zero-mean primitives.
    w = par.w;
    v    =  par.Vm*sin(w*t - th)        + par.V3*sin(3*w*t);
    dv   =  par.Vm*w*cos(w*t - th)      + 3*par.V3*w*cos(3*w*t);
    ddv  = -par.Vm*w^2*sin(w*t - th)    - 9*par.V3*w^2*sin(3*w*t);
    phi  = -par.Vm/w*cos(w*t - th)      - par.V3/(3*w)*cos(3*w*t);
    phi2 = -par.Vm/w^2*sin(w*t - th)    - par.V3/(9*w^2)*sin(3*w*t);
end

function s = phase_a(t, par)
% Hysteretic inductor in parallel with C_a. The flux linkage is imposed by
% the voltage, and the current follows the piecewise-linear loop.
    [v, dv, phi] = supply(t, par.theta(1), par);
    iL = hysteresis_current(phi, v, par);
    iC = par.Ca*dv;
    s.v = v;  s.i = iL + iC;  s.iL = iL;  s.iC = iC;  s.phi = phi;
    s.pR = zeros(size(t));                  % no resistor in this phase
    s.qL = v.*iL;                           % includes the hysteresis loss
    s.qC = v.*iC;
    s.p  = v.*s.i;
end

function i = hysteresis_current(phi, dphi, par)
% Piecewise-linear loop: slope Lsteep between the knees and Lflat elsewhere.
% Ascending branch for d(phi)/dt >= 0, descending branch otherwise.
    iA = (par.phi2 - par.phi1)/(2*par.Lflat) + (par.phi2 + par.phi1)/(2*par.Lsteep);
    iB = (par.phi2 - par.phi1)/(2*par.Lflat) - (par.phi2 + par.phi1)/(2*par.Lsteep);
    up = @(x) (x <= -par.phi1) .* (-iA + (x + par.phi2)/par.Lflat) ...
            + (x > -par.phi1 & x <= par.phi2) .* (iB + (x + par.phi1)/par.Lsteep) ...
            + (x > par.phi2) .* (iA + (x - par.phi2)/par.Lflat);
    i = zeros(size(phi));
    asc = dphi >= 0;
    i(asc)  =  up(phi(asc));
    i(~asc) = -up(-phi(~asc));
end

function s = phase_b(t, par)
% Linear parallel R_b, L_b, C_b branch in periodic steady state.
    [v, dv, phi, ddv, phi2] = supply(t, par.theta(2), par);
    s.v = v;  s.dv = dv;  s.phi = phi;
    s.iR = v/par.Rb;  s.iL = phi/par.Lb;  s.iC = par.Cb*dv;
    s.i  = s.iR + s.iL + s.iC;
    s.di = dv/par.Rb + v/par.Lb + par.Cb*ddv;           % derivative of i
    s.q  = phi/par.Rb + phi2/par.Lb + par.Cb*v;         % zero-mean primitive of i
    s.pR = v.*s.iR;  s.qL = v.*s.iL;  s.qC = v.*s.iC;  s.p = v.*s.i;
end

function s = phase_c(t, par)
% Series R_c-L_c branch connected through a transistor from ton to toff in
% each period. The current starts from zero at turn-on. After turn-off the
% inductor current circulates through the freewheeling diode and does not
% flow at the terminals.
    v   = supply(t, par.theta(3), par);
    tau = mod(t - par.ton, par.T);                      % time since the last turn-on
    on  = tau <= (par.toff - par.ton);
    i   = forced_rl(par.ton + tau, par) - forced_rl(par.ton, par).*exp(-tau*par.Rc/par.Lc);
    i(~on) = 0;
    s.v = v;  s.i = i;  s.on = on;
    s.pR = par.Rc*i.^2;
    s.qL = i.*(v - par.Rc*i);  s.qL(~on) = 0;           % L i di/dt, with L di/dt = v - R i
    s.qC = zeros(size(t));
    s.p  = v.*i;
end

function i = forced_rl(t, par)
% Steady-state response of the series R_c-L_c branch to the phase-c voltage.
    w = par.w;  th = par.theta(3);
    Z1 = hypot(par.Rc, w*par.Lc);    a1 = atan2(w*par.Lc, par.Rc);
    Z3 = hypot(par.Rc, 3*w*par.Lc);  a3 = atan2(3*w*par.Lc, par.Rc);
    i = par.Vm/Z1*sin(w*t - th - a1) + par.V3/Z3*sin(3*w*t - a3);
end

function D = duality_phase_b(B, par)
% Series branch with time-varying parameters defined by Corollary 2:
%   R_s i^2 = G v^2,  L_s i i' = Gamma v phi,  S_s i q = C v v'.
    G = 1/par.Rb;  Gam = 1/par.Lb;  Cp = par.Cb;
    v = B.v;  i = B.i;  di = B.di;  q = B.q;
    ok = abs(i) > 1e-3*max(abs(i)) & abs(di) > 1e-3*max(abs(di)) & abs(q) > 1e-3*max(abs(q));
    Rs = G*v.^2 ./ i.^2;
    Ls = Gam*v.*B.phi ./ (i.*di);
    Ss = Cp*v.*B.dv ./ (i.*q);
    D.Rs = Rs;  D.Ls = Ls;  D.Ss = Ss;  D.ok = ok;
    D.pR_series = Rs.*i.^2;  D.qL_series = Ls.*i.*di;  D.qC_series = Ss.*i.*q;
    D.kvl    = max(abs(Rs(ok).*i(ok) + Ls(ok).*di(ok) + Ss(ok).*q(ok) - v(ok))) / max(abs(v));
    D.err_pR = max(abs(D.pR_series(ok) - B.pR(ok)));
    D.err_qL = max(abs(D.qL_series(ok) - B.qL(ok)));
    D.err_qC = max(abs(D.qC_series(ok) - B.qC(ok)));
    % Energy associated with L_s(t): running integral of L_s i i', compared
    % with the energy (1/2) Gamma phi^2 stored in the constant parallel inductor.
    dt = par.T/numel(v);
    WL_int = cumtrapz(B.qL)*dt;
    WL     = 0.5*Gam*B.phi.^2;
    D.err_WL = max(abs((WL_int - WL_int(1)) - (WL - WL(1))));
    D.WLmax  = max(WL);
    D.Rs_nonneg = all(Rs(ok) >= 0);
    D.Ls_sign   = any(Ls(ok) < 0) && any(Ls(ok) > 0);
    D.Ss_sign   = any(Ss(ok) < 0) && any(Ss(ok) > 0);
end

function cmp = compare_with_reference(refdir, par)
% Relative rms difference between the waveforms computed here and the
% curves of the figures of the manuscript.
    cmp = struct();
    names = {'a', 'b', 'c'};
    model = {@phase_a, @phase_b, @phase_c};
    for k = 1:3
        file = fullfile(refdir, sprintf('powers_phase_%s.csv', names{k}));
        if ~exist(file, 'file')
            fprintf('[4-wire] reference curves not found: %s\n', file);
            continue
        end
        M = readmatrix(file, 'NumHeaderLines', 1);     % t_ms, v, i, p_active, q_L, q_C, q, p
        s = model{k}(M(:,1)*1e-3, par);
        e.v  = relrms(s.v,  M(:,2));
        e.i  = relrms(s.i,  M(:,3));
        e.qL = relrms(s.qL, M(:,5));
        e.qC = relrms(s.qC, M(:,6));
        e.p  = relrms(s.p,  M(:,8));
        e.P_model = mean(s.p(1:end-1));  e.P_ref = mean(M(1:end-1,8), 'omitnan');
        e.i_mean_ref = mean(M(1:end-1,3), 'omitnan');
        cmp.(names{k}) = e;
        fprintf('[4-wire] phase %s vs figure curves, relative rms difference: i = %.2e, p = %.2e, q_L = %.2e; P = %.1f W (figure %.1f W)\n', ...
            names{k}, e.i, e.p, e.qL, e.P_model, e.P_ref);
    end
    if isfield(cmp, 'a')
        fprintf('[4-wire] phase a: the figure curve of i_a has a mean value of %.2f A, which the symmetric loop does not reproduce\n', ...
            cmp.a.i_mean_ref);
    end
    file = fullfile(refdir, 'duality_phase_b.csv');
    if exist(file, 'file')
        M = readmatrix(file, 'NumHeaderLines', 1);     % t, Rs, Ls, Ss, p_par, p_ser, qL_par, qL_ser, qC_par, qC_ser
        cmp.duality_ref = [max(abs(M(:,5) - M(:,6)), [], 'omitnan'), ...
                           max(abs(M(:,7) - M(:,8)), [], 'omitnan'), ...
                           max(abs(M(:,9) - M(:,10)), [], 'omitnan')];
        fprintf('[4-wire] figure curves of the duality: max |parallel - series| = %.2e W (p_R), %.2e W (q_L), %.2e W (q_C)\n', ...
            cmp.duality_ref);
        Bm = phase_b(M(:,1)*1e-3, par);
        Dm = duality_phase_b(Bm, par);
        okR = Dm.ok & isfinite(M(:,2)) & abs(M(:,2)) < 50;
        okL = Dm.ok & isfinite(M(:,3)) & abs(M(:,3)) < 0.05;
        cmp.Rs_median_reldiff = median(abs(Dm.Rs(okR) - M(okR,2)) ./ abs(M(okR,2)));
        cmp.Ls_median_reldiff = median(abs(Dm.Ls(okL) - M(okL,3)) ./ abs(M(okL,3)));
        fprintf('[4-wire] series parameters vs figure curves, median relative difference: R_s = %.2e, L_s = %.2e\n', ...
            cmp.Rs_median_reldiff, cmp.Ls_median_reldiff);
    end
end

function plot_powers(par, out_pdf)
% One column per phase; each phase is plotted from the zero crossing of its voltage.
    N = 4000;
    start = [0, par.T/3, 2*par.T/3];
    model = {@phase_a, @phase_b, @phase_c};
    names = {'a', 'b', 'c'};
    fh = figure('Position', [100 100 900 760], 'Color', 'w');
    tl = tiledlayout(fh, 4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:3
        t = start(k) + (0:N).'*par.T/N;
        s = model{k}(t, par);
        tm = t*1e3;
        nexttile(tl, k);
        plot(tm, s.pR/1e3, 'b-');  grid on;  xlim([tm(1) tm(end)]);
        title(sprintf('phase $%s$', names{k}));
        if k == 1; ylabel('$p_R$ [kW]'); end
        nexttile(tl, 3 + k);
        plot(tm, s.qL/1e3, 'r-');  hold on;  plot(tm, s.qC/1e3, 'k-');  grid on;  xlim([tm(1) tm(end)]);
        if k == 1; ylabel('$q_L$, $q_C$ [kW]'); end
        if k == 3; legend({'$q_L$', '$q_C$'}, 'Location', 'northeast'); end
        nexttile(tl, 6 + k);
        plot(tm, (s.qL + s.qC)/1e3, 'r-');  grid on;  xlim([tm(1) tm(end)]);
        if k == 1; ylabel('$q$ [kW]'); end
        nexttile(tl, 9 + k);
        plot(tm, s.p/1e3, 'b-');  hold on;
        plot(tm([1 end]), mean(s.p(1:end-1))/1e3*[1 1], 'r-');  grid on;  xlim([tm(1) tm(end)]);
        xlabel('time [ms]');
        if k == 1; ylabel('$p$, $P$ [kW]'); end
        if k == 3; legend({'$p$', '$P$'}, 'Location', 'northeast'); end
    end
    exportgraphics(fh, out_pdf, 'ContentType', 'vector');
    close(fh);
end

function plot_duality(par, out_pdf)
    N = 5000;
    t = par.T/3 + (0:N).'*par.T/N;
    B = phase_b(t, par);
    D = duality_phase_b(B, par);
    tm = t*1e3;
    fh = figure('Position', [100 100 820 640], 'Color', 'w');
    tl = tiledlayout(fh, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile(tl, 1);  plot(tm, D.Rs, 'b-');  ylim([0 70]);  grid on;  xlim([tm(1) tm(end)]);
    ylabel('$R_s$ [$\Omega$]');  title('Series parameters');
    nexttile(tl, 2);  plot(tm, B.pR/1e3, 'b-');  hold on;  plot(tm, D.pR_series/1e3, 'g--');  grid on;  xlim([tm(1) tm(end)]);
    ylabel('$p_R$ [kW]');  title('Parallel (solid) and series (dashed)');
    nexttile(tl, 3);  plot(tm, D.Ls, 'r-');  ylim([-0.02 0.02]);  grid on;  xlim([tm(1) tm(end)]);
    ylabel('$L_s$ [H]');
    nexttile(tl, 4);  plot(tm, B.qL/1e3, 'r-');  hold on;  plot(tm, D.qL_series/1e3, 'g--');  grid on;  xlim([tm(1) tm(end)]);
    ylabel('$q_L$ [kW]');
    nexttile(tl, 5);  plot(tm, D.Ss/1e3, 'k-');  ylim([-3 3]);  grid on;  xlim([tm(1) tm(end)]);
    ylabel('$S_s$ [mF$^{-1}$]');  xlabel('time [ms]');
    nexttile(tl, 6);  plot(tm, B.qC/1e3, 'k-');  hold on;  plot(tm, D.qC_series/1e3, 'g--');  grid on;  xlim([tm(1) tm(end)]);
    ylabel('$q_C$ [kW]');  xlabel('time [ms]');
    exportgraphics(fh, out_pdf, 'ContentType', 'vector');
    close(fh);
end

function e = relrms(x, ref)
    ok = isfinite(x) & isfinite(ref);
    e = sqrt(mean((x(ok) - ref(ok)).^2)) / sqrt(mean(ref(ok).^2));
end

function y = rms_local(x)
    y = sqrt(mean(x.^2));
end
