function write_manuscript_tables
% Export the three reconstruction tables, three supplementary tables, and GCV table.
% Main tables bold ties at the displayed precision and retain all stop categories.
folder=fileparts(mfilename('fullpath')); out=fullfile(folder,'output','tables');
if ~exist(out,'dir'), mkdir(out); end
families={'gaussian','defocus','motion'};
order=[1 3 2]; names={'No log-det','Proposed','Marg'}; stops={'Limit','Grad.','Step'};
for k=1:3
    a=load(fullfile(folder,'output',['reconstruction_' families{k}],'results.mat'));
    main={'\begin{tabular}{@{}lclrrrl@{}}';'\toprule'; ...
        'Method & Initial Guess & $\by^*$ & RRE & SSIM & PSF err. & Stop \\';'\midrule'};
    extra={'\begin{tabular}{@{}lcrrr@{}}';'\toprule'; ...
        'Method & \shortstack{\small Initial\\Guess} & PE & RR & PSNR \\';'\midrule'};
    best=[min(round(a.RRE(:),3)),max(round(a.SSIM(:),3)),min(round(a.PSF_error(:),3))];
    best_extra=[min(round(a.parameter_error(:),3)),max(round(a.PSNR(:),2))];
    for block=1:3
        m=order(block);
        for start=1:3
            y=a.y_estimate(:,m,start); parts=cell(numel(y),1);
            for j=1:numel(y)
                if y(j)~=0 && abs(y(j))<1e-3, parts{j}=scientific(y(j));
                else, parts{j}=sprintf('%.3f',y(j)); end
            end
            estimate=strjoin(parts,',');
            if numel(y)>1, estimate=['(' estimate ')']; end
            line=[names{m} ' & ' num2str(start) ' & $' estimate '$ & ' ...
                metric(a.RRE(m,start),best(1)) ' & ' metric(a.SSIM(m,start),best(2)) ...
                ' & ' metric(a.PSF_error(m,start),best(3)) ' & ' stops{a.stop(m,start)+1} ' \\'];
            main{end+1,1}=line; %#ok<AGROW>
            extra{end+1,1}=[names{m} ' & ' num2str(start) ' & $' ...
                extra_metric(a.parameter_error(m,start),best_extra(1),3) '$ & $' ...
                scientific(a.relative_residual(m,start)) '$ & $' ...
                extra_metric(a.PSNR(m,start),best_extra(2),2) '$ \\']; %#ok<AGROW>
        end
        if block<3
            main{end+1,1}='\addlinespace[3pt]'; extra{end+1,1}='\addlinespace[3pt]'; %#ok<AGROW>
        end
    end
    main=[main;{'\bottomrule';'\end{tabular}'}];
    extra=[extra;{'\bottomrule';'\end{tabular}'}];
    main_text=strjoin(main,newline); extra_text=strjoin(extra,newline);
    write_text(fullfile(out,['tabular_reconstruction-' families{k} '.tex']),main_text);
    write_text(fullfile(out,['tabular_additional-' families{k} '.tex']),extra_text);
end
a=load(fullfile(folder,'output','gcv_single_sample','results.mat'));
lines={'\begin{tabular}{lrrrrrr}';'\toprule'; ...
    'Method & $\alpha$ & $\lambda$ & $\by^*$ & RRE & PSF error & SSIM\\';'\midrule'};
names={'No log-det','Automatic selection'};
for j=1:2
    lines{end+1,1}=sprintf('%s & %.3f & $%.3g$ & $(%.3f,%.3f)$ & %.3f & %.3f & %.3f \\\\', ...
        names{j},a.alpha(j),a.lambda(j),a.y_estimate(:,j),a.RRE(j),a.PSF_error(j),a.SSIM(j)); %#ok<AGROW>
end
lines=[lines;{'\bottomrule';'\end{tabular}'}]; tabular=strjoin(lines,newline);
write_text(fullfile(out,'tabular_auto.tex'),tabular);
fprintf('Exported the seven numerical tables to %s\n',out);

end

function s=scientific(value)
parts=strsplit(sprintf('%.2e',value),'e');
s=[parts{1} '\times10^{' num2str(str2double(parts{2})) '}'];
end
function s=metric(value,best)
s=sprintf('%.3f',value);
if round(value,3)==best, s=['\mathbf{' s '}']; end
s=['$' s '$'];
end
function s=extra_metric(value,best,precision)
s=sprintf('%.*f',precision,value);
if round(value,precision)==best, s=['\mathbf{' s '}']; end
end
function write_text(file,text)
fid=fopen(file,'w'); assert(fid>=0,'Cannot write %s',file);
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',text);
end
