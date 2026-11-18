function paradox_resistive_unbalanced()
% PARADOX_RESISTIVE_UNBALANCED  §IV.B (rewritten, May 2026)
%   Single-branch resistive three-wire load (de Leon-Cohen scenario, 2005):
%   only R_ab between phases a and b; phases b-c and c-a open.
%   Demonstrates that
%     (i) branch-level decomposition: p_term = p_R,ab; dW_L = dW_C = 0;
%     (ii) Czarnecki CPC active currents are nonzero in all 3 phases,
%          including phase c whose physical current is identically zero;
%     (iii) Akagi p-q reports nonzero q_alphabeta for a purely resistive load.

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    p2_figstyle();
    % Local style overrides (Fig 4 only): heavier lines + larger labels
    set(groot, 'defaultAxesFontSize', 20);
    set(groot, 'defaultTextFontSize', 20);
    set(groot, 'defaultLegendFontSize', 21);
    set(groot, 'defaultLineLineWidth', 2.6);
    figdir = fullfile(here, '..', 'results', 'figures');
    resdir = fullfile(here, '..', 'results');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    if ~exist(resdir,'dir'); mkdir(resdir); end

    f0 = 50;
    om = 2*pi*f0;
    fs = 20e3;
    t  = (0:1/fs:3/f0).';
    Vm = 400*sqrt(2);

    % Balanced source line-to-line voltages
    vab = Vm*cos(om*t);
    vbc = Vm*cos(om*t - 2*pi/3);
    vca = -(vab + vbc);

    % Single resistor between phases a and b (rest of delta open)
    R_ab = 0.65;
    iab = vab/R_ab;
    ibc = zeros(size(t));
    ica = zeros(size(t));

    ia =  iab - ica;     % = iab
    ib =  ibc - iab;     % = -iab
    ic =  ica - ibc;     % = 0

    % Virtual neutral
    van = (vab - vca)/3;
    vbn = (vbc - vab)/3;
    vcn = (vca - vbc)/3;

    p_term = van.*ia + vbn.*ib + vcn.*ic;
    p_R    = R_ab*iab.^2;
    fprintf('[Resistive open] residual p_term - p_R: %.3e (relative: %.3e)\n', ...
        norm(p_term - p_R), norm(p_term - p_R)/norm(p_term));

    % Akagi p-q
    K = sqrt(2/3);
    valpha = K*(van - 0.5*vbn - 0.5*vcn);
    vbeta  = K*(0 + sqrt(3)/2*vbn - sqrt(3)/2*vcn);
    ialpha = K*(ia - 0.5*ib - 0.5*ic);
    ibeta  = K*(0 + sqrt(3)/2*ib - sqrt(3)/2*ic);
    p_pq = valpha.*ialpha + vbeta.*ibeta;
    q_pq = valpha.*ibeta - vbeta.*ialpha;

    % CPC active currents (using equivalent conductance G_e)
    % G_e = P / (||v||^2_3phase),  where P is mean active power and
    % ||v||^2 = mean(v_an^2 + v_bn^2 + v_cn^2).
    P_avg = mean(p_term);
    Vnorm2 = mean(van.^2 + vbn.^2 + vcn.^2);
    G_e = P_avg / Vnorm2;

    ia_cpc_act = G_e*van;
    ib_cpc_act = G_e*vbn;
    ic_cpc_act = G_e*vcn;   % nonzero even though physical i_c = 0

    fprintf('[Resistive open] G_e = %.4f S; CPC active RMS per phase: %.3f, %.3f, %.3f A\n', ...
        G_e, rms(ia_cpc_act), rms(ib_cpc_act), rms(ic_cpc_act));
    fprintf('[Resistive open] physical phase c current RMS: %.3e A (zero by construction)\n', rms(ic));
    fprintf('[Resistive open] q_pq RMS = %.2f var (nonzero on a purely resistive load)\n', rms(q_pq));

    fh = figure('Position',[100 100 680 720],'Color','w');

    subplot(4,1,1);
    co = lines(3);
    yyaxis left; plot(t*1e3, vab,'-','Color',co(1,:)); hold on; plot(t*1e3, vbc,'-','Color',co(2,:)); plot(t*1e3, vca,'-','Color',co(3,:));
    ylabel('voltages [V]');
    yyaxis right; plot(t*1e3, ia,'--','Color',co(1,:)); plot(t*1e3, ib,'--','Color',co(2,:)); plot(t*1e3, ic,'--','Color',co(3,:));
    ylabel('currents [A]');
    ax=gca; ax.YAxis(1).Color='k'; ax.YAxis(2).Color='k';
    legend({'$v_{ab}$','$v_{bc}$','$v_{ca}$','$i_a$','$i_b$','$i_c\equiv 0$'}, ...
        'Location','northeast','NumColumns',2);
    title('Single-branch resistive load: $R_{ab}=0.65\,\Omega$, phase $c$ open');
    grid on;

    subplot(4,1,2);
    plot(t*1e3, p_term/1e3,'-'); hold on;
    plot(t*1e3, p_R/1e3,'--');
    legend({'$p_{\rm term}$','$p_{R,ab}$'},'Location','northeast');
    ylabel('[kW]');
    title('Branch-level decomposition: $p_{\rm term}=p_{R,ab}$, $W_L^{\prime}=W_C^{\prime}=0$');
    grid on;

    subplot(4,1,3);
    plot(t*1e3, ia_cpc_act,'-'); hold on;
    plot(t*1e3, ib_cpc_act,'-');
    plot(t*1e3, ic_cpc_act,'-','LineWidth',2);
    plot(t*1e3, ic,':k','LineWidth',1.5);
    legend({'$i_{a,{\rm CPC,act}}$','$i_{b,{\rm CPC,act}}$', ...
            '$i_{c,{\rm CPC,act}}$ (phantom)','$i_c$ physical $\equiv 0$'}, ...
        'Location','northeast','NumColumns',2);
    ylabel('[A]');
    title('CPC active currents populate every phase, including the open one');
    grid on;

    subplot(4,1,4);
    plot(t*1e3, p_pq/1e3); hold on; plot(t*1e3, q_pq/1e3);
    legend({'$p_{\alpha\beta}$ [kW]','$q_{\alpha\beta}$ [kvar]'},'Location','northeast');
    xlabel('time [ms]'); ylabel('Akagi $p$--$q$');
    title('Akagi $q_{\alpha\beta}\neq 0$ for a purely resistive load');
    grid on;

    out_pdf = fullfile(figdir,'fig_p2_resistive_unbalanced_pq.pdf');
    set(findall(fh,'Type','Legend'),'FontSize',21);   % enforce larger legend text
    exportgraphics(fh, out_pdf,'ContentType','vector');
    close(fh);

    save(fullfile(resdir,'paradox_resistive_unbalanced.mat'), ...
        't','vab','vbc','vca','ia','ib','ic','iab','p_term','p_R', ...
        'p_pq','q_pq','ia_cpc_act','ib_cpc_act','ic_cpc_act','R_ab','G_e');

    fprintf('[Resistive open] figure: %s\n', out_pdf);
end
