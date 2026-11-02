function verify_tables
% Compare all seven generated table bodies with the paper reference tables.
folder=fileparts(mfilename('fullpath'));
reference=dir(fullfile(folder,'reference_tables','*.tex'));
assert(numel(reference)==7,'Expected seven manuscript reference tables.');
for k=1:numel(reference)
    expected=fileread(fullfile(reference(k).folder,reference(k).name));
    actual=fileread(fullfile(folder,'output','tables',reference(k).name));
    expected=regexprep(expected,'\s',''); actual=regexprep(actual,'\s','');
    assert(strcmp(actual,expected),'Table differs from the paper reference: %s',reference(k).name);
end
fprintf('PASS: all seven tables match the paper references (ignoring whitespace).\n');
end
