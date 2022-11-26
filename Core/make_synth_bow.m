% pathstr = '/Users/mikhail/Desktop/LOCAL_DATA/pcacoreset/';
% prefix = 'synth_bow_';

num_blocks = 1;

%%
for block_no = 1:num_blocks

    %filename = [prefix num2str(block_no)];
    %fprintf('block %d: %s\n',block_no,filename);
    
    ni = n/num_blocks;
    Ai = SVDCoresetTest.gen_noisy_subspace(ni,d,r,j);
    x1 = (block_no-1)*ni+1;
    x2 = block_no*ni;
    
    A = sparse(n,d);
    A(x1:x2,:) = Ai;
    
    try
        [~] = A(n,d);
    catch
        A(n,d) = 0;
    end 
    
    %save('-v7.3',[pathstr output_subdir prefix num2str(block_no)],'A','block_no')
    %save([pathstr prefix num2str(block_no)],'A','block_no')
    
end
