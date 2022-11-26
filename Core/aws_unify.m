% cluster1 = parcluster('mikhail_c3x64_1');
% cluster2 = parcluster('mikhail_c3x64_2');

k = 10;

%%
load(sprintf('data/k%d/block1',k))
unified_coreset = S.D;
unified_coreset_error = unified_coreset.coreset_errors(end);
cumulative_errors = [];
cumulative_sizes = [];
num_points_streamed = [];
cumulative_errors(1) = unified_coreset_error;
cumulative_sizes(1) = unified_coreset.coreset_size;
num_points_streamed(1) = S.stream.num_points_streamed;
processed_blocks = zeros(1,64);
processed_blocks(S.block_no) = 1;
save(sprintf('data/k%d/wiki_unified_k%d',k,k), ...
    'processed_blocks','num_points_streamed','unified_coreset','unified_coreset_error','cumulative_errors','cumulative_sizes')

%%
for i = 2:64
    
    load(sprintf('data/k%d/wiki_unified_k%d',k,k))
    
    if not(processed_blocks(i))
        i
        load(sprintf('data/k%d/block%d',k,i))
        unified_coreset = S.stream.coreset_alg.mergedCoreset(unified_coreset.C,S.D.C);
        unified_coreset_error = unified_coreset.coreset_errors(end);
        cumulative_errors(S.block_no) = unified_coreset_error;
        cumulative_sizes(S.block_no) = unified_coreset.coreset_size;
        processed_blocks(S.block_no) = 1;
        num_points_streamed = [num_points_streamed num_points_streamed(end)+S.stream.num_points_streamed];
        assert(length(num_points_streamed)==i)
        save(sprintf('data/k%d/wiki_unified_k%d',k,k), ...
            'processed_blocks','num_points_streamed','unified_coreset','unified_coreset_error','cumulative_errors','cumulative_sizes')
    end
    
    if not(processed_blocks(i))
        break
    end
    
end

