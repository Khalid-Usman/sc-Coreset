setenv('PATH','/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:/usr/local/sbin')

format compact
format short g

project_paths = {pwd};

disp('Adding project paths:')
if ispc
    error('Platform not supported!')
elseif ismac
    project_paths = cat(1,project_paths,'/Users/mikhail/Documents/MIT/pcacoreset');
    project_paths = cat(1,project_paths,'/Users/mikhail/MIT/DATA/pcacoreset');
    project_paths = cat(1,project_paths,'/Users/mikhail/Desktop/LOCAL_DATA/pcacoreset/');
elseif isunix && not(ismac)
    error('Platform not supported!')
end
project_paths

for i = 1:length(project_paths)
    paths_str = genpath(project_paths{i});
    paths = strsplit(paths_str,':')';
    used = false(1,length(paths));
    for j = 1:length(paths)
        svn_dir = not(isempty(strfind(paths{j},'.svn')));
        unused_dir = not(isempty(strfind(paths{j},'__unused__')));
        used(j) = not(svn_dir) & not(unused_dir);
    end
    paths = strjoin(paths(used),':');
    addpath(paths)
end
disp('Done!')

clear project_paths i

myclear
