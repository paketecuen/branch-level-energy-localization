function paradox_topology_indeterminacy()
% PARADOX_TOPOLOGY_INDETERMINACY  Section IV.D.
%   Constructive proof of Theorem 2.
%   Same terminal data (v_ab, v_bc, i_a, i_b) is reproduced by two
%   admissible three-wire equivalents:
%     (M1) delta-parallel RL with branch parameters (R_xy, L_xy)
%     (M2) wye-series virtual-neutral RL with branch parameters (R_x, L_x)
%   Both reconstruct p_term within numerical precision; the branch
%   energy decompositions differ.

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    p2_figstyle();
    set(groot, 'defaultAxesFontSize', 18);
    set(groot, 'defaultTextFontSize', 18);
    set(groot, 'defaultLegendFontSize', 17);
    set(groot, 'defaultLineLineWidth', 2.4);
    figdir = fullfile(here, '..', 'results', 'figures');
    resdir = fullfile(here, '..', 'results');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    if ~exist(resdir,'dir'); mkdir(resdir); end

    f0 = 50;
    om = 2*pi*f0;
    fs = 50e3;
    Tend = 3/f0;
    t  = (0:1/fs:Tend - 1/fs).';   % integer-cycle window for exact averaging

    Vm = 400*sqrt(2);

    % --- Truth (M1): unbalanced delta-parallel RL ---
    Rab = 4.0;  Lab = 0.010;
    Rbc = 6.0;  Lbc = 0.012;
    Rca = 5.0;  Lca = 0.011;

    % Pure fundamental balanced source
    vab = Vm*cos(om*t);
    vbc = Vm*cos(om*t - 2*pi/3);
    vca = -(vab + vbc);
    vab_d = -Vm*om*sin(om*t);
    vbc_d = -Vm*om*sin(om*t - 2*pi/3);
    vca_d = -(vab_d + vbc_d);

    % Branch admittances (parallel RL): Y(om) = 1/R + 1/(j om L)
    Yab = 1/Rab - 1j/(om*Lab);
    Ybc = 1/Rbc - 1j/(om*Lbc);
    Yca = 1/Rca - 1j/(om*Lca);

    % Line-to-line phasors
    Vab_ph = Vm;
    Vbc_ph = Vm*exp(-1j*2*pi/3);
    Vca_ph = -(Vab_ph + Vbc_ph);

    % Branch currents (phasor, then time)
    Iab_ph = Yab*Vab_ph;
    Ibc_ph = Ybc*Vbc_ph;
    Ica_ph = Yca*Vca_ph;

    iab = real(Iab_ph*exp(1j*om*t));
    ibc = real(Ibc_ph*exp(1j*om*t));
    ica = real(Ica_ph*exp(1j*om*t));

    ia = iab - ica;
    ib = ibc - iab;
    ic = ica - ibc;

    % Branch RL currents (split conductive vs inductive)
    iab_R = vab/Rab;       iab_L = iab - iab_R;
    ibc_R = vbc/Rbc;       ibc_L = ibc - ibc_R;
    ica_R = vca/Rca;       ica_L = ica - ica_R;

    % Analytical derivatives of inductor currents (for energy term)
    iab_L_d = (1/Lab)*vab;     % since L*di/dt = v on inductor branch (parallel)
    ibc_L_d = (1/Lbc)*vbc;
    ica_L_d = (1/Lca)*vca;

    % Branch energy terms (delta description)
    pR_ab = Rab*iab_R.^2;     % equivalent to vab^2/Rab
    pR_bc = Rbc*ibc_R.^2;
    pR_ca = Rca*ica_R.^2;
    pR_delta = pR_ab + pR_bc + pR_ca;

    Wdot_ab = Lab*iab_L.*iab_L_d;
    Wdot_bc = Lbc*ibc_L.*ibc_L_d;
    Wdot_ca = Lca*ica_L.*ica_L_d;
    Wdot_delta = Wdot_ab + Wdot_bc + Wdot_ca;

    % Virtual neutral
    van = (vab - vca)/3;
    vbn = (vbc - vab)/3;
    vcn = (vca - vbc)/3;
    p_term = van.*ia + vbn.*ib + vcn.*ic;
    p_recon_delta = pR_delta + Wdot_delta;
    eps_delta = norm(p_term - p_recon_delta)/norm(p_term);
    fprintf('[Topology indeterminacy] delta reconstruction residual: %.3e\n', eps_delta);

    % --- M2: wye-virtual series RL fit using phasors at fundamental ---
    % v_xn'(t) = R_x i_x(t) + L_x di_x/dt
    % At fundamental: V_xn = (R_x + j om L_x) * I_x  =>  Z_x = V_xn / I_x

    Van_ph = (Vab_ph - Vca_ph)/3;
    Vbn_ph = (Vbc_ph - Vab_ph)/3;
    Vcn_ph = (Vca_ph - Vbc_ph)/3;

    Ia_ph = Iab_ph - Ica_ph;
    Ib_ph = Ibc_ph - Iab_ph;
    Ic_ph = Ica_ph - Ibc_ph;

    Za = Van_ph/Ia_ph;  Ra = real(Za);  La = imag(Za)/om;
    Zb = Vbn_ph/Ib_ph;  Rb = real(Zb);  Lb = imag(Zb)/om;
    Zc = Vcn_ph/Ic_ph;  Rc = real(Zc);  Lc = imag(Zc)/om;

    fprintf('[Topology indeterminacy] Wye-virtual: Ra=%.3f Rb=%.3f Rc=%.3f Ohm\n', Ra,Rb,Rc);
    fprintf('[Topology indeterminacy]              La=%.3f Lb=%.3f Lc=%.3f mH\n', La*1e3,Lb*1e3,Lc*1e3);

    % Time-domain reconstruction with M2
    iad = -Ia_ph*om;       % di/dt at fundamental: factor j*om in phasor; time = real(j*om*Ia*exp(...))
    ia_d = real(1j*om*Ia_ph*exp(1j*om*t));
    ib_d = real(1j*om*Ib_ph*exp(1j*om*t));
    ic_d = real(1j*om*Ic_ph*exp(1j*om*t));

    pR_a = Ra*ia.^2;  pR_b = Rb*ib.^2;  pR_c = Rc*ic.^2;
    pR_wye = pR_a + pR_b + pR_c;
    Wdot_a = La*ia.*ia_d;
    Wdot_b = Lb*ib.*ib_d;
    Wdot_c = Lc*ic.*ic_d;
    Wdot_wye = Wdot_a + Wdot_b + Wdot_c;
    p_recon_wye = pR_wye + Wdot_wye;
    eps_wye = norm(p_term - p_recon_wye)/norm(p_term);
    fprintf('[Topology indeterminacy] wye-virtual reconstruction residual: %.3e\n', eps_wye);

    % Sanity: difference between recovered branch energy maps
    fprintf('[Topology indeterminacy] mean(pR_delta)=%.1f W; mean(pR_wye)=%.1f W\n', ...
        mean(pR_delta), mean(pR_wye));
    fprintf('[Topology indeterminacy] mean(W_delta) profile range %.1f W; mean(W_wye) range %.1f W\n', ...
        max(Wdot_delta)-min(Wdot_delta), max(Wdot_wye)-min(Wdot_wye));

    % --- Windowed LS identification of branch parameters ---
    % Delta (parallel RL, branch xy):  i_xy = G_xy * v_xy + Gamma_xy * ∫v_xy dt
    % Wye-virtual (series RL, phase x): v_xn' = R_x * i_x + L_x * di_x/dt
    Twin_id = 1/f0;   % one fundamental period
    IntVab = cumtrapz(t, vab);
    IntVbc = cumtrapz(t, vbc);
    IntVca = cumtrapz(t, vca);
    [G_ab_t, Gamma_ab_t] = wlsq2(iab, vab, IntVab, t, Twin_id);
    [G_bc_t, Gamma_bc_t] = wlsq2(ibc, vbc, IntVbc, t, Twin_id);
    [G_ca_t, Gamma_ca_t] = wlsq2(ica, vca, IntVca, t, Twin_id);

    dia_dt = gradient(ia, t);
    dib_dt = gradient(ib, t);
    dic_dt = gradient(ic, t);
    [Ra_t, La_t] = wlsq2(van, ia, dia_dt, t, Twin_id);
    [Rb_t, Lb_t] = wlsq2(vbn, ib, dib_dt, t, Twin_id);
    [Rc_t, Lc_t] = wlsq2(vcn, ic, dic_dt, t, Twin_id);

    fprintf('[Topology indeterminacy] Windowed identification (median of valid window):\n');
    fprintf('  Delta: G_ab=%.4f, G_bc=%.4f, G_ca=%.4f S (truth %.4f, %.4f, %.4f)\n', ...
        median(G_ab_t,'omitnan'), median(G_bc_t,'omitnan'), median(G_ca_t,'omitnan'), ...
        1/Rab, 1/Rbc, 1/Rca);
    fprintf('  Delta: Gamma_ab=%.2f, Gamma_bc=%.2f, Gamma_ca=%.2f 1/H (truth %.2f, %.2f, %.2f)\n', ...
        median(Gamma_ab_t,'omitnan'), median(Gamma_bc_t,'omitnan'), median(Gamma_ca_t,'omitnan'), ...
        1/Lab, 1/Lbc, 1/Lca);
    fprintf('  Wye:   R_a=%.3f, R_b=%.3f, R_c=%.3f Ohm\n', ...
        median(Ra_t,'omitnan'), median(Rb_t,'omitnan'), median(Rc_t,'omitnan'));
    fprintf('  Wye:   L_a=%.3f, L_b=%.3f, L_c=%.3f mH\n', ...
        median(La_t,'omitnan')*1e3, median(Lb_t,'omitnan')*1e3, median(Lc_t,'omitnan')*1e3);

    % Figure
    fh = figure('Position',[100 100 760 760],'Color','w');
    tl = tiledlayout(fh,4,2,'TileSpacing','tight','Padding','tight');

    nexttile(tl);
    plot(t*1e3, vab); hold on; plot(t*1e3, vbc); plot(t*1e3, vca);
    legend({'$v_{ab}$','$v_{bc}$','$v_{ca}$'},'Location','northeast','NumColumns',3);
    ylabel('voltage [V]');
    title('Line-to-line voltages');
    grid on;

    nexttile(tl);
    plot(t*1e3, ia); hold on; plot(t*1e3, ib); plot(t*1e3, ic);
    legend({'$i_a$','$i_b$','$i_c$'},'Location','northeast','NumColumns',3);
    ylabel('current [A]');
    title('Line currents');
    grid on;

    nexttile(tl);
    plot(t*1e3, pR_ab/1e3); hold on;
    plot(t*1e3, pR_bc/1e3);
    plot(t*1e3, pR_ca/1e3);
    legend({sprintf('$p_{R,ab}$, $R{=}%.1f\\,\\Omega$', Rab), ...
            sprintf('$p_{R,bc}$, $R{=}%.1f\\,\\Omega$', Rbc), ...
            sprintf('$p_{R,ca}$, $R{=}%.1f\\,\\Omega$', Rca)}, ...
        'Location','northeast');
    ylabel('Joule [kW]');
    title('M1: delta dissipation');
    grid on;

    nexttile(tl);
    plot(t*1e3, pR_a/1e3); hold on;
    plot(t*1e3, pR_b/1e3);
    plot(t*1e3, pR_c/1e3);
    legend({sprintf('$p_{R,a}$, $R{=}%.2f\\,\\Omega$', Ra), ...
            sprintf('$p_{R,b}$, $R{=}%.2f\\,\\Omega$', Rb), ...
            sprintf('$p_{R,c}$, $R{=}%.2f\\,\\Omega$', Rc)}, ...
        'Location','northeast');
    ylabel('Joule [kW]');
    title('M2: wye-virtual dissipation');
    grid on;

    nexttile(tl);
    plot(t*1e3, Wdot_ab/1e3); hold on;
    plot(t*1e3, Wdot_bc/1e3);
    plot(t*1e3, Wdot_ca/1e3);
    legend({sprintf('$W_{L,ab}^{\\prime}$, $L{=}%.0f$ mH', Lab*1e3), ...
            sprintf('$W_{L,bc}^{\\prime}$, $L{=}%.0f$ mH', Lbc*1e3), ...
            sprintf('$W_{L,ca}^{\\prime}$, $L{=}%.0f$ mH', Lca*1e3)}, ...
        'Location','northeast');
    ylabel('storage rate [kW]');
    title('M1: delta storage rate');
    grid on;

    nexttile(tl);
    plot(t*1e3, Wdot_a/1e3); hold on;
    plot(t*1e3, Wdot_b/1e3);
    plot(t*1e3, Wdot_c/1e3);
    legend({sprintf('$W_{L,a}^{\\prime}$, $L{=}%.2f$ mH', La*1e3), ...
            sprintf('$W_{L,b}^{\\prime}$, $L{=}%.2f$ mH', Lb*1e3), ...
            sprintf('$W_{L,c}^{\\prime}$, $L{=}%.2f$ mH', Lc*1e3)}, ...
        'Location','northeast');
    ylabel('storage rate [kW]');
    title('M2: wye-virtual storage rate');
    grid on;

    % --- Row 4: identified parameter time evolution (windowed LS) ---
    co = lines(3);
    nexttile(tl);
    yyaxis left;
    plot(t*1e3, G_ab_t,'-', 'Color', co(1,:), 'LineWidth', 2.4); hold on;
    plot(t*1e3, G_bc_t,'-', 'Color', co(2,:), 'LineWidth', 2.4);
    plot(t*1e3, G_ca_t,'-', 'Color', co(3,:), 'LineWidth', 2.4);
    ylabel('$G_{xy}(t)$ [S]');
    set(gca,'YColor','k');
    ylim([0.13 0.29]);
    yyaxis right;
    plot(t*1e3, Gamma_ab_t,'--', 'Color', co(1,:), 'LineWidth', 2.4); hold on;
    plot(t*1e3, Gamma_bc_t,'--', 'Color', co(2,:), 'LineWidth', 2.4);
    plot(t*1e3, Gamma_ca_t,'--', 'Color', co(3,:), 'LineWidth', 2.4);
    set(gca,'YColor','k');
    ylim([75 110]);
    xlim([10 50]); xlabel('time [ms]');
    legend({'$G_{ab}$','$G_{bc}$','$G_{ca}$', ...
            '$\Gamma_{ab}$','$\Gamma_{bc}$','$\Gamma_{ca}$'}, ...
        'Location','northeast','NumColumns',2);
    title('M1: delta parameters');
    grid on;

    nexttile(tl);
    yyaxis left;
    plot(t*1e3, Ra_t,'-', 'Color', co(1,:), 'LineWidth', 2.4); hold on;
    plot(t*1e3, Rb_t,'-', 'Color', co(2,:), 'LineWidth', 2.4);
    plot(t*1e3, Rc_t,'-', 'Color', co(3,:), 'LineWidth', 2.4);
    set(gca,'YColor','k');
    ylim([0.40 0.65]);
    yyaxis right;
    plot(t*1e3, La_t*1e3,'--', 'Color', co(1,:), 'LineWidth', 2.4); hold on;
    plot(t*1e3, Lb_t*1e3,'--', 'Color', co(2,:), 'LineWidth', 2.4);
    plot(t*1e3, Lc_t*1e3,'--', 'Color', co(3,:), 'LineWidth', 2.4);
    ylabel('$L_x(t)$ [mH]');
    set(gca,'YColor','k');
    ylim([2.0 2.85]);
    xlim([10 50]); xlabel('time [ms]');
    legend({'$R_a$','$R_b$','$R_c$','$L_a$','$L_b$','$L_c$'}, ...
        'Location','northeast','NumColumns',2);
    title('M2: wye-virtual parameters');
    grid on;

    out_pdf = fullfile(figdir,'fig_p2_threewire_indeterminacy.pdf');
    set(findall(fh,'Type','Legend'),'FontSize',18);   % enforce larger legend text
    exportgraphics(fh, out_pdf,'ContentType','vector');
    close(fh);

    save(fullfile(resdir,'paradox_topology_indeterminacy.mat'), ...
        't','vab','vbc','vca','ia','ib','ic','p_term', ...
        'pR_ab','pR_bc','pR_ca','Wdot_ab','Wdot_bc','Wdot_ca', ...
        'pR_a','pR_b','pR_c','Wdot_a','Wdot_b','Wdot_c', ...
        'Ra','Rb','Rc','La','Lb','Lc','Rab','Rbc','Rca','Lab','Lbc','Lca');

    fprintf('[Topology indeterminacy] figure: %s\n', out_pdf);
end

% ----------------------------------------------------- windowed LS helper
function [p1_t, p2_t] = wlsq2(target, reg1, reg2, t, Twin)
% Sliding-window least-squares fit  target(idx) ≈ p1*reg1(idx) + p2*reg2(idx).
%   Window of duration Twin centered on each interior sample. Within-window
%   means are removed from each column to absorb DC offsets (e.g. of cumulative
%   integrals) so that p1, p2 are identified consistently.
    fs = 1/(t(2)-t(1));
    Nw = round(Twin*fs);
    if mod(Nw,2)==0, Nw = Nw+1; end
    half = (Nw-1)/2;
    N = length(t);
    p1_t = nan(size(t));
    p2_t = nan(size(t));
    for k = (half+1):(N-half)
        idx = (k-half):(k+half);
        r1 = reg1(idx) - mean(reg1(idx));
        r2 = reg2(idx) - mean(reg2(idx));
        b  = target(idx) - mean(target(idx));
        x  = [r1, r2] \ b;
        p1_t(k) = x(1);
        p2_t(k) = x(2);
    end
end
