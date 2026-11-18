function paradox_czarnecki_CL()
% PARADOX_CZARNECKI_CL  Single-phase pedagogical illustration of Theorem 2.
%   Replicates Czarnecki (1997) parallel L=32/57 H, C=49/32 F under a
%   non-sinusoidal voltage. Shows that branch-level energy localization
%   distinguishes inductive from capacitive storage even when fundamental
%   currents alone do not.
%
% Outputs:
%   figures_paper2/fig_p2_czarnecki_CL_illustration.pdf
%   results/paradox_czarnecki_CL.mat

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    p2_figstyle();
    figdir = fullfile(here, '..', 'results', 'figures');
    resdir = fullfile(here, '..', 'results');
    if ~exist(figdir,'dir'); mkdir(figdir); end
    if ~exist(resdir,'dir'); mkdir(resdir); end

    % Parameters from Czarnecki 1997
    f  = 1/(2*pi);
    om = 2*pi*f;          % base angular frequency (=1 rad/s)
    C  = 49/32;
    L  = 32/57;

    % Time grid: two periods of fundamental
    T  = 1/f;
    N  = 4000;
    t  = linspace(0, 2*T, N).';

    % Excitation: v(t) = sqrt(2)*(100 sin om t + 25 sin 3 om t)
    v   = sqrt(2)*( 100*sin(om*t) + 25*sin(3*om*t) );
    vd  = sqrt(2)*( 100*om*cos(om*t) + 25*3*om*cos(3*om*t) );

    % Capacitor branch: i_C = C*dv/dt
    iC  = C*vd;

    % Inductor branch (analytic): i_L = (1/L) integral v dt -- with zero mean
    iL  = sqrt(2)*( 100/(om*L)*sin(om*t - pi/2) + 25/(3*om*L)*sin(3*om*t - pi/2) );
    iLd = sqrt(2)*( 100/(om*L)*om*cos(om*t - pi/2) + 25/(3*om*L)*3*om*cos(3*om*t - pi/2) );

    i   = iC + iL;

    % Energy terms
    p       = v.*i;
    Wdot_L  = L*iL.*iLd;
    Wdot_C  = C*v.*vd;
    p_recon = Wdot_L + Wdot_C;
    residual = p - p_recon;

    eps_e = norm(residual)/norm(p);
    fprintf('[Czarnecki C/L] |residual|/|p| = %.3e\n', eps_e);

    fh = figure('Position',[100 100 680 620],'Color','w');

    subplot(3,1,1);
    yyaxis left;  plot(t, v); ylabel('$v$ [V]');
    yyaxis right; plot(t, i); ylabel('$i$ [A]');
    title('Terminal waveforms');
    grid on; xlim([0 2*T]);

    subplot(3,1,2);
    plot(t, p,'-'); hold on;
    plot(t, p_recon,'--');
    legend({'$p(t)=v\,i$','$W_L^{\prime}+W_C^{\prime}$'}, 'Location','best');
    ylabel('power [W]');
    title(sprintf('Terminal power and branch reconstruction ($\\varepsilon_e=%.1e$)', eps_e));
    grid on; xlim([0 2*T]);

    subplot(3,1,3);
    plot(t, Wdot_L); hold on;
    plot(t, Wdot_C);
    legend({'$W_L^{\prime}=L\,i_L i_L^{\prime}$','$W_C^{\prime}=C\,v v^{\prime}$'}, 'Location','best');
    ylabel('storage rate [W]'); xlabel('time [s]');
    title('Branch storage rates: distinct profiles');
    grid on; xlim([0 2*T]);

    out_pdf = fullfile(figdir, 'fig_p2_czarnecki_CL_illustration.pdf');
    exportgraphics(fh, out_pdf, 'ContentType','vector');
    close(fh);

    save(fullfile(resdir,'paradox_czarnecki_CL.mat'), ...
        't','v','i','iC','iL','p','Wdot_L','Wdot_C','residual','eps_e','C','L','om');

    fprintf('[Czarnecki C/L] figure: %s\n', out_pdf);
end
