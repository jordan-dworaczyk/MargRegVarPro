function run_all(mode)
% RUN_ALL Recompute the paper, or export figures from locally generated results.
% run_all / run_all('full'): all experiments, weight searches, figures, tables.
% run_all('experiments'): skip the three reference-image weight searches.
% run_all('figures'): regenerate figures/tables from existing output/*.mat.
if nargin<1, mode='full'; end
mode=validatestring(mode,{'full','experiments','figures'});
folder=fileparts(mfilename('fullpath')); previous_path=path;
cleanup=onCleanup(@() path(previous_path)); %#ok<NASGU>
addpath(folder,'-begin');
verify_code;
stages={'run_objective_landscapes'};
if ~strcmp(mode,'figures')
    stages=[stages,{'run_reconstruction_gaussian','run_reconstruction_defocus', ...
        'run_reconstruction_motion','run_auxiliary_scale_sensitivity','run_gcv_single_sample'}];
    if strcmp(mode,'full')
        stages=[stages,{'run_alpha_sweep_gaussian','run_alpha_sweep_defocus','run_alpha_sweep_motion'}];
    end
end
stages=[stages,{'plot_ground_truth_reference','plot_blurred_data', ...
    'plot_reconstruction_composites','plot_reconstruction_landscapes', ...
    'plot_auxiliary_scale_sensitivity','plot_gcv_single_sample','write_manuscript_tables'}];
for k=1:numel(stages)
    fprintf('\n[%d/%d] %s\n',k,numel(stages),stages{k});
    started=tic;
    feval(stages{k});
    fprintf('Completed %s in %.1f seconds.\n',stages{k},toc(started));
end
verify_tables;
fprintf('PASS: %s workflow completed.\n',mode);
end
