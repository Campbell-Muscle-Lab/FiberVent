function create_single_summary_file(options)
% Create a summary file showing system properties at end
% of simulation

arguments
    options.top_data_folder = ""
    options.sim_file_name = "sim_output.txt";

    options.last_n_points = 2500;
    options.movmean_window = 20
    options.rel_peak_prominence = 0.75
    options.t_nan_s = 0.15

    options.extract_mode = "diastolic";
    options.vent_chamber_radius_field = "vent_chamber_radius";
    options.keep_fields = [ ...
        "time", ...
        "vent_wall_thickness", "vent_wall_volume", ...
        "vent_chamber_radius", "vent_chamber_height", ...
        "vent_ejection_fraction", ...
        "vent_efficiency", "vent_ATP_used_per_s", ...
        "vent_stroke_work_J", "vent_stroke_energy_used_J", ...
        "mito_volume", ...
        "mito_ATP_generated_M_per_liter_per_s", ...
        "hr_heart_rate_bpm"];

    options.summary_file_string = ""

    options.figure_zoom = 1;
    options.figure_zoom_folder = "";

    options.setup_file_string
    options.experiment_mode = "";
    options.multiplier_variables = [];
end

% Code

% Find data files in top folder
files = dir(fullfile(options.top_data_folder, "**", ...
    options.sim_file_name));
data_files = string(fullfile({files.folder}, {files.name}))';

% Sort them in natural order
if (numel(data_files) > 1)
    data_files = sort_nat(data_files);
end

% Loop through them, saving data to out as we go
out = [];
w = waitbar(0);

for i = 1 : numel(data_files)

    % Progress
    waitbar(i/numel(data_files), w, ...
        sprintf('File %i of %i', i, numel(data_files)));

    % Read in a file
    d = readtable(data_files(i));
    dn = d.Properties.VariableNames';

    % Find the last n summary points
    sd = size(d);
    vi = find(d.summary_full_dump == 1);
    vi = vi(numel(vi)-options.last_n_points:end);

    d = d(vi, :);
    dn = d.Properties.VariableNames';

    % Save some data
    % Get the file number
    [file_path, file_name] = fileparts(data_files(i));
    out.file_path(i) = file_path;
    out.file_name(i) = file_name;

    file_parts = strsplit(data_files(i), filesep);
    out.condition_number(i) = str2double(file_parts(end-1));

    r = d.(options.vent_chamber_radius_field);
    if (options.extract_mode == "diastolic")
        [out.vent_chamber_radius(i), ki] = max(r);
    else
        [out.vent_chamber_radius(i), ki] = min(r);
    end
    for j = 1 : numel(options.keep_fields)
        out.(options.keep_fields(j))(i) = d.(options.keep_fields(j))(ki);
    end

    % Pull off some other stuff
    if (options.figure_zoom > 0)

        dt = max(diff(d.time));
        d.pressure_0_smooth = smoothdata(d.pressure_0, ...
            "movmean", options.movmean_window);
        d.dp0dt_smooth = [0 ; diff(d.pressure_0_smooth) ./ dt];
        d.hsl_smooth = smoothdata(d.fs_hs_length, ...
            "movmean", options.movmean_window);
        d.dhsldt_smooth = [0 ; diff(d.hsl_smooth) ./ dt];

        % Find peaks in the dpdt
        dpdt_stats = summary_stats(d.dp0dt_smooth);
        [dpdt_max, dpdt_max_ind] = findpeaks(d.dp0dt_smooth, ...
            MinPeakProminence = ...
            options.rel_peak_prominence * dpdt_stats.range);

        % Look for troughs in the dpdt, which is tricky because
        % we have to watch out for transients
        dpdt_min = [];
        dpdt_min_ind = [];
        if (numel(dpdt_max_ind) > 0)
            dpdt_temp = d.dp0dt_smooth;
            % Nan before the first max
            dpdt_temp(1:dpdt_max_ind(1)) = NaN;
            % NaN a portion after the max
            for j = 1 : numel(dpdt_max_ind)
                ti_nan = find( (d.time >= d.time(dpdt_max_ind(j))) & ...
                                (d.time <= (d.time(dpdt_max_ind(j)) + ...
                                    options.t_nan_s)) );
                dpdt_temp(ti_nan) = NaN;
            end
            % Nan after the last max
            dpdt_temp(dpdt_max_ind(end):end) = NaN;
            % Find the troughs, again hard
            nan_ind = isnan(dpdt_temp);
            vi = find(diff(nan_ind));
            t_vi = d.time(vi);
            for j = 1 : 2 : numel(t_vi)
                y_temp = dpdt_temp(vi(j):vi(j+1));
                [a,b] = min(y_temp);
                dpdt_min = [dpdt_min a];
                dpdt_min_ind = [dpdt_min_ind vi(j) + b - 1];
            end
        end

        % Find peaks and troughs in the hsl
        hsl_stats = summary_stats(d.hsl_smooth);
        [hsl_max, hsl_max_ind] = findpeaks(d.hsl_smooth, ...
            MinPeakProminence = ...
            options.rel_peak_prominence * hsl_stats.range);

        [hsl_min, hsl_min_ind] = findpeaks(-d.hsl_smooth, ...
            MinPeakProminence = ...
            options.rel_peak_prominence * hsl_stats.range);
        hsl_min = -hsl_min;

        % Find peaks and troughs in the dhsl
        dhsl_stats = summary_stats(d.dhsldt_smooth);
        [dhsldt_max, dhsldt_max_ind] = findpeaks(d.dhsldt_smooth, ...
            MinPeakProminence = ...
            options.rel_peak_prominence * dhsl_stats.range);

        [dhsldt_min, dhsldt_min_ind] = findpeaks(-d.dhsldt_smooth, ...
            MinPeakProminence = ...
            options.rel_peak_prominence * dhsl_stats.range);
        dhsldt_min = -dhsldt_min;

        % Save data

        % First set defaults
        out.dpdt_max(i) = NaN;
        out.dpdt_min(i) = NaN;
        out.hsl_max(i) = NaN;
        out.hsl_min(i) = NaN;
        out.dhsldt_max(i) = NaN;
        out.dhsldt_min(i) = NaN;

        % Now fill them in
        if (numel(dpdt_max) > 0)
            out.dpdt_max(i) = mean(dpdt_max);
        end

        if (numel(dpdt_min) > 0)
            out.dpdt_min(i) = mean(dpdt_min);
        end

        if (numel(hsl_max) > 0)
            out.hsl_max(i) = mean(hsl_max);
        end

        if (numel(hsl_min) > 0)
            out.hsl_min(i) = mean(hsl_min);
        end

        if (numel(dhsldt_max) > 0)
            out.dhsldt_max(i) = mean(dhsldt_max);
        end

        if (numel(dhsldt_min) > 0)
            out.dhsldt_min(i) = mean(dhsldt_min);
        end

        % Make the figure
        [sp, fig_options] = ...
            layout_subplots( ...
                figure_handle = options.figure_zoom, ...
                figure_width = 7, ...
                panels_wide = 2, ...
                panels_high = 2, ...
                padding_left = 0.75, ...
                padding_right = 0.25, ...
                x_to_y_ratio = 2);

        subplot(sp(1));
        hold on;
        plot(d.time, d.pressure_0, 'b-');
        plot(d.time, d.pressure_0_smooth, 'r-');

        subplot(sp(3));
        hold on;
        plot(d.time, d.dp0dt_smooth, 'r-');
        if (~isnan(out.dpdt_max(i)))
            plot(d.time([1 end]), out.dpdt_max(i) * [1 1], 'g--');
            subplot(sp(1));
            plot(d.time(dpdt_max_ind), d.pressure_0_smooth(dpdt_max_ind), ...
                'go');
        end
        if (~isnan(out.dpdt_max(i)))
            subplot(sp(3))
            plot(d.time([1 end]), out.dpdt_min(i) * [1 1], 'g--');
            subplot(sp(1));
            plot(d.time(dpdt_min_ind), d.pressure_0_smooth(dpdt_min_ind), ...
                'go');
        end

        subplot(sp(2));
        hold on;
        plot(d.time, d.fs_hs_length, 'b-');
        plot(d.time, d.hsl_smooth, 'r-');
        if (~isnan(out.hsl_max(i)))
            plot(d.time([1 end]), out.hsl_max(i) * [1 1], 'g--');
        end
        if (~isnan(out.hsl_min(i)))
            plot(d.time([1 end]), out.hsl_min(i) * [1 1], 'g--');
        end

        subplot(sp(4));
        hold on;
        plot(d.time, d.dhsldt_smooth, 'r-');
        if (~isnan(out.dhsldt_max(i)))
            plot(d.time([1 end]), out.dhsldt_max(i) * [1 1], 'g--');
            subplot(sp(2));
            plot(d.time(dhsldt_max_ind), d.hsl_smooth(dhsldt_max_ind), ...
                'go');
        end
        if (~isnan(out.dhsldt_min(i)))
            subplot(sp(4))
            plot(d.time([1 end]), out.dhsldt_min(i) * [1 1], 'g--');
            subplot(sp(2));
            plot(d.time(dhsldt_min_ind), d.hsl_smooth(dhsldt_min_ind), ...
                'go');
        end

        % Save figure
        if (~isempty(options.figure_zoom_folder))
            % Prep the folder
            if (~exist('figure_zoom_created'))
                if (isdir(options.figure_zoom_folder))
                    % Clean it
                    delete(fullfile(options.figure_zoom_folder, '*'));
                else
                    % Make it
                    mkdir(options.figure_zoom_folder);
                end
                figure_zoom_created = 1;
            end

            % Now write the figure to disk
            output_file_string = fullfile(options.figure_zoom_folder, ...
                sprintf('%s_cond_%i.png', ...
                    out.file_name(i), ...
                    out.condition_number(i)));

            disp(sprintf('Saving to: %s', output_file_string));

            exportgraphics(fig_options.figure_handle, ...
                output_file_string, ...
                Resolution = 600);
        end
    end

end

% Convert data to table format
out = columnize_structure(out);
out = struct2table(out);

% Now we need to add information about the peturbation
s = readstruct(options.setup_file_string);

if (options.experiment_mode == "multipliers")

    omv = options.multiplier_variables{1}
    n_m = numel(omv)

    for i = 1 : n_m
        % Find the right adjustment
        adj = s.FiberVent_setup.model.manipulations.adjustments;
        for j = 1 : numel(adj)

            if (adj(j).variable == omv(i))
                adj_ind = j;
                break
            end
        end

        % Now loop through the summary table
        class_string = sprintf('class_%i', i);
        variable_string = sprintf('variable_%i', i);
        multiplier_string = sprintf('multiplier_%i', i);

        for j = 1 : size(out, 1)
            out.(class_string)(j) = s.FiberVent_setup.model.manipulations.adjustments(adj_ind).class;
            out.(variable_string)(j) = s.FiberVent_setup.model.manipulations.adjustments(adj_ind).variable;
            out.(multiplier_string)(j) = s.FiberVent_setup.model.manipulations.adjustments(adj_ind). ...
                multipliers(out.condition_number(j));
        end
    end

elseif (options.experiment_mode == "perturbations")

    % Scan through the table
    for i = 1 : size(out, 1)
        cond_ind = out.condition_number(i);
            
        % Find the reight perturbations
        pert = s.FiberVent_setup.characterization(1).perturbation;
        pert_ind = [];
        for j = 1 : numel(pert)
            if (pert(j).simulation(1) == cond_ind)
                pert_ind = [pert_ind j];
            end
        end

        % More than one
        for j = 1 : numel(pert_ind)
            temp_class = sprintf('class_%i', j);
            temp_variable = sprintf('variable_%i', j);
            temp_total_change = sprintf('total_change_%i', j);

            out.(temp_class)(i) = pert(pert_ind(j)).class;
            out.(temp_variable)(i) = pert(pert_ind(j)).variable;
            out.(temp_total_change)(i) = pert(pert_ind(j)).total_change;
        end
    end
end

% Move variables
out_n = out.Properties.VariableNames;
move_ind = find( strcmp(out_n, 'condition_number') | ...
                    startsWith(out_n, 'class') | ...
                    startsWith(out_n, 'variable') | ...
                    startsWith(out_n, 'total_change') | ...
                    startsWith(out_n, "multiplier") );

out = movevars(out, out_n(move_ind), "Before", 1)

% Write table to file
try
    delete(options.summary_file_string)
end
writetable(out, options.summary_file_string)

