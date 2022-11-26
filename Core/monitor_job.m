function monitor_job(cluster,job_ID,timeout)

while(1)
    clc
    job = cluster.findJob('ID',job_ID);
    diary(job)
    pause(timeout)
end

% function monitor_job(jobs,timeout)
% if not(iscell(jobs))
%     jobs = {jobs};
% end
% if nargin < 2
%     timeout = 1;
% end
% while(1)
%     for i = 1:length(jobs)
%         clc
%         job = jobs{i};
%         if not(isempty(job))
%             fprintf('----- Job %d: ID = %d -----\n\n',i,job.ID)
%             diary(job)
%             fprintf('\n----- Job %d: ID = %d -----\n',i,job.ID)
%             fprintf('State = %s\n',job.State)
%             fprintf('Start Time = %s\n',job.StartTime)
%             pause(timeout)
%         end
%     end
% end
