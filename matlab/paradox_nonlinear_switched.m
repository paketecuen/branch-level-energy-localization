function paradox_nonlinear_switched()
% PARADOX_NONLINEAR_SWITCHED  Section IV.C
%   Three-wire delta-connected TRIAC-switched resistive load.
%   Each delta branch (a-b, b-c, c-a) has a TRIAC + R in series and fires
%   alpha after the zero-crossing of its own line-to-line voltage; it
%   conducts until the next zero-crossing (resistive load: i = v/R).
%   The line currents are differences of switching branch currents and
%   show heavy distortion. Classical theories (Akagi p-q, Fryze) report
%   large nonactive components although the load contains no storage.
%   Branch-level identification returns G_xy(t) toggling between 0 and
%   1/R, with Gamma_xy = 0 and C_xy = 0 everywhere.

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    p2_figstyle();
    set(groot, 'defaultAxesFontSize', 20);
    set(groot, 'defaultTextFontSize', 20);
    set(groot, 'defaultLegendFontSize', 21);
    set(groot, 'defaultLineLineWidth', 2.6);
    figdir = fullfile(here, '..', 'results', 'figures');
    resdir = fullfile(here, '..', 'results');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    if ~exist(resdir,'dir'); mkdir(resdir); end

    f0 = 50; om = 2*pi*f0; fs = 50e3;
    Tend = 3/f0;
    t  = (0:1/fs:Tend).';
    Vm_LL = 400*sqrt(2);  % line-to-line peak (V_LL,rms = 400 V)
    R     = 10;
    alpha = pi/3;

    % Balanced line-to-line voltages
    vab = Vm_LL*cos(om*t);
    vbc = Vm_LL*cos(om*t - 2*pi/3);
    vca = Vm_LL*cos(om*t + 2*pi/3);

    % Virtual neutral (centroid)
    van = (vab - vca)/3;
    vbn = (vbc - vab)/3;
    vcn = (vca - vbc)/3;

    % Each delta branch: TRIAC + R in series, fires on its own line-line voltage
    sab = triac_switch_local(vab, t, om, alpha);
    sbc = triac_switch_local(vbc, t, om, alpha);
    sca = triac_switch_local(vca, t, om, alpha);
    iab = sab.*vab/R;
    ibc = sbc.*vbc/R;
    ica = sca.*vca/R;

    % Line currents
    ia = iab - ica;
    ib = ibc - iab;
    ic = ica - ibc;

    % Branch energy: only Joule, no storage
    pR_ab = R*iab.^2; pR_bc = R*ibc.^2; pR_ca = R*ica.^2;
    pR_tot = pR_ab + pR_bc + pR_ca;
    p_term = van.*ia + vbn.*ib + vcn.*ic;

    % Classical scalar indicators
    Vrms_phase = sqrt(mean(van.^2 + vbn.^2 + vcn.^2));
    Irms_line  = sqrt(mean(ia.^2 + ib.^2 + ic.^2));
    P    = mean(p_term);
    S_aggr = Vrms_phase*Irms_line;
    QF   = sqrt(max(S_aggr^2 - P^2, 0));

    % Akagi p-q from phase quantities
    K = sqrt(2/3);
    valpha = K*(van - 0.5*vbn - 0.5*vcn);
    vbeta  = K*(0 + sqrt(3)/2*vbn - sqrt(3)/2*vcn);
    ialpha = K*(ia - 0.5*ib - 0.5*ic);
    ibeta  = K*(0 + sqrt(3)/2*ib - sqrt(3)/2*ic);
    p_pq = valpha.*ialpha + vbeta.*ibeta;
    q_pq = valpha.*ibeta - vbeta.*ialpha;

    % Identified branch conductance (parallel R-only topology)
    G_ab = iab ./ vab;  G_ab(abs(vab)<1) = NaN;
    G_bc = ibc ./ vbc;  G_bc(abs(vbc)<1) = NaN;
    G_ca = ica ./ vca;  G_ca(abs(vca)<1) = NaN;

    fprintf('[Nonlinear switched, 3-wire delta] P=%.1f W, S=%.1f VA, Q_F=%.1f VAr\n', ...
        P, S_aggr, QF);
    fprintf('[Nonlinear switched] residual (p_term - sum p_R) rel: %.3e\n', ...
        norm(p_term - pR_tot)/norm(p_term));

    fh = figure('Position',[100 100 680 720],'Color','w');

    subplot(4,1,1);
    yyaxis left;  plot(t*1e3, vab); ylabel('$v_{ab}$ [V]');
    yyaxis right; plot(t*1e3, iab,'-'); hold on;
    plot(t*1e3, ia,':');
    ylabel('current [A]');
    legend({'$v_{ab}$','$i_{ab}$ branch','$i_a$ line'}, ...
        'Location','northeast');
    title(sprintf('3-wire delta TRIAC: branch $a$-$b$ and line $a$ ($\\alpha=%.0f^\\circ$)', ...
        alpha*180/pi));
    grid on;

    subplot(4,1,2);
    plot(t*1e3, p_pq/1e3); hold on;
    plot(t*1e3, q_pq/1e3);
    legend({'$p_{\alpha\beta}$ [kW]','$q_{\alpha\beta}$ [kvar]'}, ...
        'Location','northeast');
    ylabel('Akagi $p$--$q$');
    title(sprintf('Classical $q_{\\alpha\\beta}\\neq 0$; aggregate $Q_F=%.0f$~VAr', QF));
    grid on;

    subplot(4,1,3);
    plot(t*1e3, p_term/1e3,'-'); hold on;
    plot(t*1e3, pR_tot/1e3,'--');
    legend({'$p_{\rm term}$','$\sum p_R$'}, ...
        'Location','northeast');
    ylabel('[kW]');
    title('Terminal power equals total Joule dissipation (no storage)');
    grid on;

    subplot(4,1,4);
    plot(t*1e3, pR_ab/1e3); hold on;
    plot(t*1e3, pR_bc/1e3);
    plot(t*1e3, pR_ca/1e3);
    legend({'$p_{R,ab}$','$p_{R,bc}$','$p_{R,ca}$'}, ...
        'Location','northeast');
    xlabel('time [ms]'); ylabel('per-branch [kW]');
    ymax = ceil(max([pR_ab; pR_bc; pR_ca])/1e3/5)*5 + 5;
    ylim([0 ymax]);
    title('Branch-localized Joule dissipation ($\Gamma_{xy}\equiv 0,\,C_{xy}\equiv 0$)');
    grid on;

    out_pdf = fullfile(figdir,'fig_p2_nonlinear_switched.pdf');
    set(findall(fh,'Type','Legend'),'FontSize',21);   % enforce larger legend text
    exportgraphics(fh, out_pdf,'ContentType','vector');
    close(fh);

    save(fullfile(resdir,'paradox_nonlinear_switched.mat'), ...
        't','vab','vbc','vca','van','vbn','vcn','ia','ib','ic', ...
        'iab','ibc','ica','p_term','pR_tot','p_pq','q_pq', ...
        'G_ab','G_bc','G_ca','P','S_aggr','QF','R','alpha');
    fprintf('[Nonlinear switched] figure: %s\n', out_pdf);
end

% ------------------------------------------------------------- TRIAC mask
function s = triac_switch_local(v, t, om, alpha)
% TRIAC conduction mask for a sinusoidal voltage v(t).
%   Recovers the sin-form phase phi at t(1) from v(1) and dv/dt(1) so that
%   v ≈ A*sin(phi(t)). The TRIAC fires alpha after each zero-crossing of v
%   and conducts until the next zero-crossing (resistive load: current
%   crosses zero with v). Mask = 1 when mod(phi, pi) > alpha.
    A = max(abs(v));
    if A == 0 || numel(v) < 2
        s = zeros(size(t));
        return
    end
    dt = t(2) - t(1);
    v0    = v(1);
    dvdt0 = (v(2) - v(1))/dt;
    phi0  = atan2(v0/A, dvdt0/(A*om));
    phi_t = phi0 + om*(t - t(1));
    s = double(mod(phi_t, pi) > alpha);
end
