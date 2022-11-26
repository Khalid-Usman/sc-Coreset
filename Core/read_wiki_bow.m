pathstr = '~/MIT/DATA/pcacoreset/wiki_en1/'

input_subdir = 'bow_mm_blocks/';
output_subdir = 'bow_mat/';
prefix = 'wiki_en_bow_';

files = dir([pathstr input_subdir]);
isdir = cell2mat(extractfield(files,'isdir'));
filenames = extractfield(files,'name');
filenames = filenames(not(isdir));

header_fid = fopen([pathstr input_subdir filenames{1}]);
[~] = fgets(header_fid);
tline = fgets(header_fid);
sz = str2num(tline);
n = sz(1)
d = sz(2)
s = sz(3)/(n*d)

filenames = filenames(2:end);

%%
for block_no = 1:length(filenames)
    
%     i = floor((block_no-1)/26+1);
%     j = mod(block_no-1,26)+1;
%     suffix = char([96+i 96+j]);
%     filename = ['bow_mm_blocks/' prefix suffix];

    filename = filenames{block_no};
    fprintf('block %d: %s\n',block_no,filename);
    
    filepath = [pathstr input_subdir filename];
    [is,js,ss] = import_mm_block(filepath);
    A = sparse(is,js,ss);
    
    try
        [~] = A(sz(1),sz(2));
    catch
        A(sz(1),sz(2)) = 0;
    end 
    
    save([pathstr output_subdir prefix num2str(block_no)],'A','block_no')
    
end
