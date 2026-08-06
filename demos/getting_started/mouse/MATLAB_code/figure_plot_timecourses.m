function figure_plot_timecourses(options)
% Plot simulation time-courses

arguments
    options.data_files = [ ...
        "../sim_data/sim_output/1/sim_output.txt", ...
        "../sim_data/sim_output/7/sim_output.txt"]
    
    options.template_files = [
        "../templates/template_example_beats.json"]
    
    options.output_image_files = [
        "../output/fig_fibervent_mouse"]
    
    options.envelope_no_of_bins = 50

    options.output_image_types = ["png", "svg"]

    options.data_filter = struct("column", "time", "type", "isfinite");
end

% Set the colormap
cm = return_matplotlib_default_colors();
cm = repmat(cm, [1 1 2]);
cm(:, :, 2) = brighten(cm(:, :, 1), 0.5)

% Make the figure
figure_multi_x( ...
    options.data_files, ...
    options.template_files(1), ...
    data_file_filters = options.data_filter, ...
    output_file_string = options.output_image_files(1), ...
    output_file_types = options.output_image_types, ...
    envelope_no_of_bins = options.envelope_no_of_bins, ...
    trace_color_map = cm);
