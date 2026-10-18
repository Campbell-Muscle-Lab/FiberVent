function figure_fibervent_adaptation(options)

    arguments
        options.sim_results_file_string = "../sim_data/sim_output/1/sim_output.txt"

        options.status_folder = "../sim_data/sim_output/1/status";

    end

    if (1)
    
        % Get the status files
        status_files = sort_nat(findfiles('json', options.status_folder))';
    
        % Initialize the profile
        z = [];
    
        waitbar(0);
        % Load the iso-type as a function of length
        for fi = 1 : numel(status_files)
    
            waitbar(fi/numel(status_files));
    
            y = return_isotype_profile(status_files{fi}, "m");
            z = [z y'];
        end

        save('temp.mat');
    else
        a = load("temp.mat");

        z = a.z;
    end


    % Average
    n  = 10;
    for ci = 1 : size(z, 2)/n
        vi = ((ci-1) * n) + (1:n);
        z2(:, ci) = mean(z(:, vi), 2);
    end
    z = z2;

    z = smoothdata(z, "movmean", 10);

    figure(1);
    clf
    waterfall(z')

    