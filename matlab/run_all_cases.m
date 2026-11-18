function run_all_cases()
% RUN_ALL_CASES  Regenerates the waveforms, balance residuals and figures of
%   every test case of the manuscript. Outputs are written to ../results
%   (figures in ../results/figures, console log in ../results/run_log.txt).

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    resdir = fullfile(here, '..', 'results');
    if ~exist(resdir, 'dir'); mkdir(resdir); end
    logfile = fullfile(resdir, 'run_log.txt');
    if exist(logfile, 'file'); delete(logfile); end
    wstate = warning('off', 'MATLAB:print:ContentTypeImageSuggested');
    tbstate = get(groot, 'defaultAxesToolbarVisible');
    set(groot, 'defaultAxesToolbarVisible', 'off');     % keep the axes toolbar out of the exported figures
    diary(logfile);

    fprintf('Branch-level energy localization in three-phase loads: test cases\n');
    fprintf('MATLAB %s\n\n', version);

    paradox_czarnecki_CL();            % Section II-G, single-phase illustration
    paradox_balanced_3ph_RL();         % Section IV-A
    paradox_balanced_RL_compare();     % Section IV-A, comparison with IEEE Std. 1459
    paradox_resistive_unbalanced();    % Section IV-B
    paradox_nonlinear_switched();      % Section IV-C
    paradox_topology_indeterminacy();  % Section IV-D
    paradox_fluctuating_phase();       % Section IV-E
    fourwire_nonlinear();              % Section IV-F

    fprintf('\nDone.\n');
    diary('off');
    warning(wstate);
    set(groot, 'defaultAxesToolbarVisible', tbstate);

    % The scripts print absolute paths of the figures; keep them relative in the log.
    txt = fileread(logfile);
    txt = strrep(txt, [here filesep '..' filesep], '');
    fid = fopen(logfile, 'w');  fwrite(fid, txt);  fclose(fid);
end
