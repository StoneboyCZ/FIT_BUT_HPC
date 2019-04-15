function simulateRun(nRLCSizes, workersNum, runsNum)
%% simulateRun 
% runs the simulateTime script with Telegraph Line Problem 
% PARAMETERS:
% nRLCSizes - array with numbers of segments (e.g. nRLCSizes =
% 50:100:1000;)
% runsNum - number of runs for each number of segments
% workersNum - number of runs for each number of segments
%
% NOTE: MTSM_one_serial and MTSM_one_parallel also includes the
% precalculation time (see script simulateTime.m, lines 229 - 231)
%
% NOTE2: simulateTime.m - close all and clc was deleted to keep printouts

    all_MTSM_serial = 0.0;
    all_MTSM_one_serial = 0.0;
    all_MTSM_one_parallel = 0.0;

    for i=1:length(workersNum)
        fprintf('\n*** TELEGRAPH LINE PROBLEM : runs = %d, workers = %d *** \n', runsNum, workersNum(i));
        for j=1:length(nRLCSizes)
            fprintf('nRLC: %d\n',nRLCSizes(j))
            
            for k=1:runsNum
                [MTSM_serial, MTSM_one_serial, MTSM_one_parallel] = simulateExportTimes(nRLCSizes(j),workersNum(i));
                all_MTSM_serial = all_MTSM_serial + MTSM_serial;
                all_MTSM_one_serial = all_MTSM_one_serial + MTSM_one_serial;
                all_MTSM_one_parallel = all_MTSM_one_parallel + MTSM_one_parallel;      
            end
            
            fprintf('Average time: Serial MTSM solver: %d seconds \n', all_MTSM_serial/runsNum);
            fprintf('Average time: Serial MTSM solver - one mtx: %d seconds \n', all_MTSM_one_serial/runsNum);
            fprintf('Average time: Parallel MTSM solver - one mtx: %d seconds \n', all_MTSM_one_parallel/runsNum);
            fprintf('---------------------------\n');

            all_MTSM_serial = 0.0;
            all_MTSM_one_serial = 0.0;
            all_MTSM_one_parallel = 0.0;
        end
    end
    
    
%     for i=1:length(nRLCSizes)
%         fprintf('\n*** TELEGRAPH LINE PROBLEM : nRLC = %d, runs = %d *** \n', nRLCSizes(i), runsNum);
%         for j=1:runsNum
%             for k=1:length(workersNum)
%                 fprintf('%d workers\n',workersNum(k));
%                 parpool('local',workersNum(k))
%                 [MTSM_serial, MTSM_one_serial, MTSM_one_parallel] = simulateExportTimes(nRLCSizes(i));
%                 all_MTSM_serial = all_MTSM_serial + MTSM_serial;
%                 all_MTSM_one_serial = all_MTSM_one_serial + MTSM_one_serial;
%                 all_MTSM_one_parallel = all_MTSM_one_parallel + MTSM_one_parallel;
%             end
%         end
%         
%         fprintf('Average time: Serial MTSM solver: %d seconds \n', all_MTSM_serial/runsNum);
%         fprintf('Average time: Serial MTSM solver - one mtx: %d seconds \n', all_MTSM_one_serial/runsNum);
%         fprintf('Average time: Parallel MTSM solver - one mtx: %d seconds \n', all_MTSM_one_parallel/runsNum);
%         fprintf('---------------------------\n');
% 
%         all_MTSM_serial = 0.0;
%         all_MTSM_one_serial = 0.0;
%         all_MTSM_one_parallel = 0.0;
%     end

    % shut down the parallel pool
    poolobj = gcp('nocreate');
    delete(poolobj);
end