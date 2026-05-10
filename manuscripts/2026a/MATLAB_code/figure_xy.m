function figure_xy(options)

    arguments

        options.table_folder = "../summaries"

        options.sarcomere_fields = ["membrames", "thin_parameters", ...
            "titin_parameters", "myofilament_kinetics"];

        options.max_heart_rate = 150
        options.scan_parameters_min_time = 640

        options.marker_size = 25;
        options.marker_face_alpha = 0.25;
        options.marker_edge_color = 0.5 *  ones(1,3)
        options.marker_edge_alpha = 0.5;

        options.marker_mitochondria = 's'
        options.marker_size_mitochondria = 25
        options.marker_face_color_mitochondria = [1 0 0]
        options.marker_face_alpha_mitochondria = 1
        options.marker_edge_color_mitochondria = [0 0 0]

        options.marker_contractility = 'o'
        options.marker_size_contractility = 25
        options.marker_face_color_contractility = [1 1 1];
        options.marker_face_alpha_contractility = 0.5
        options.marker_edge_color_contractility = 0.7 * ones(1,3)

        options.marker_hypertension = '^'
        options.marker_size_hypertension = 25
        options.marker_face_color_hypertension = [0 0.75 0];
        options.marker_face_alpha_hypertension = 1
        options.marker_edge_color_hypertension = 0 * ones(1,3)

        options.marker_av_leak = 'v'
        options.marker_size_av_leak = 25
        options.marker_face_color_av_leak = [0 0 1];
        options.marker_face_alpha_av_leak = 1
        options.marker_edge_color_av_leak = 0 * ones(1,3)

        options.marker_mv_leak_clean = 'd'
        options.marker_size_mv_leak_clean = 25
        options.marker_face_color_mv_leak_clean = [1 0.6 0];
        options.marker_face_alpha_mv_leak_clean = 1
        options.marker_edge_color_mv_leak_clean = 0 * ones(1,3)

        options.marker_mv_leak_reduced = 'd'
        options.marker_size_mv_leak_reduced = 15
        options.marker_face_color_mv_leak_reduced = [0 0 0];
        options.marker_face_alpha_mv_leak_reduced = 0.5
        options.marker_edge_color_mv_leak_reduced = 0 * ones(1,3)

        options.annotation_font_size = 10

        options.base_annotation_x = 33
        options.base_annotation_y = 9
        options.base_annotation_text_x_offset = 3.5;
        options.base_annotation_text_y_offset = -0.6;
        options.base_annotation_label = "Base" + newline + " model"
        options.base_annotation_color = 0 * ones(1,3);

        options.mybpc_annotation = true
        options.mybpc_annotation_x = 30
        options.mybpc_annotation_y = 12
        options.mybpc_annotation_color = 0 * ones(1,3);
        options.mybpc_annotation_text_x_offset = 1.7;
        options.mybpc_annotation_text_y_offset = 0.25;
        options.mybpc_annotation_label = ...
            "Less stabilization" + newline + "of myosin OFF" + ...
            newline + "by MyBP-C"

        options.lab_1 = [ ...
            ["Myofilament /" + newline + "Ca^{2+}-handling" + ...
                newline + "variants" + newline + "colored by" + ...
                newline + "dPdt"] ; ...
            ["Mitochondrial" + newline + "dysfunction"] ]

        options.lab_2 = [ ...
            ["Hypertension"] ; ...
            ["AV insufficency"] ; ...
            ["MV insufficiency"] ; ...
            ["MV insufficiency" + newline + "and depressed" + newline + ...
                "contractility" + newline + "colored by" + ...
                newline + "dPdt"] ];

        options.legend_font_size = 8
        options.legend_icon_column_width = 10
        options.legend_position = [1.05 0.94]
        options.legend_vertical_scale = 1.2

        options.legend_position_2 = [1.05 0.91]
        options.legend_vertical_scale_2 = 1.3

        options.x_scaling_factor = 1000;
        options.y_scaling_factor = 1000;
        
        options.x_ticks = [20 40]
        options.y_ticks = [5 14]

        options.tick_font_size = 11

        options.x_label_offset = -0.2;
        options.y_label_offset = -0.32;

        options.c_limits = [1000 4000];
        options.c_limits_2 = [1000 3000];

        options.title_strings = { ...
            {'Growth associated with', 'altered myocyte properties'}, ...
            {'Growth associated with', 'hemodynamic perturbations'} };
        options.title_y_offset = 1.2

        options.output_image_file = "../figures/fig_xy"
        options.output_image_types = ["png", "svg"];
    end

    % Load data
    all_data_files = findfiles('xlsx', options.table_folder)'
    d = readtable(all_data_files{1});
    dn = d.Properties.VariableNames;

    if (numel(all_data_files) > 1)
        for i = 2 : numel(all_data_files)
            % Read in new table and make sure columns match
            dn = d.Properties.VariableNames;
            d2 = readtable(all_data_files{i});
            d2n = d2.Properties.VariableNames;
            ndiff = [setdiff(dn, d2n) setdiff(d2n, dn)];

            for j = 1 : numel(ndiff)
                if (~any(strcmp(dn, ndiff{j})))
                    if (iscell(d2.(ndiff{j})))
                        d.(ndiff{j}) = repmat({''}, [size(d, 1) 1]);
                    else
                        d.(ndiff{j}) = NaN * ones(size(d,1), 1);
                    end
                    n = numel(d.Properties.VariableNames);
                else
                    if (iscell(d.(ndiff{j})))
                        d2.(ndiff{j}) = repmat({''}, [size(d2, 1), 1]);
                    else
                        d2.(ndiff{j}) = NaN * ones(size(d2,1), 1);
                    end
                    n2 = numel(d2.Properties.VariableNames);
                end
            end
            
            % Now that they are matching, collate
            d = [d ; d2];
        end
    end

    % Make figure
    [subplots, fig_options] = layout_subplots( ...
        figure_width = 7, ...
        panels_wide = 2, ...
        x_to_y_ratio = 0.65, ...
        padding_left = 0.8, ...
        padding_right = 1.25, ...
        padding_top = 0.75, ...
        padding_bottom = 0.75);

    % XY
    subplot(subplots(1));

    % Find the base model for later annotation
    vi = find( (d.total_change_1 == 0) & ...
                (d.total_change_2 == 0) & ...
                (d.total_change_3 == 0) );

    base_radius = d.vent_chamber_radius(vi);
    base_thickness = d.vent_wall_thickness(vi);

    % Find and plot perturbations that are sarcomere fields alone
    vi = [];
    for i = 1 : size(d,1)
        if (any(startsWith(d.class_1(i), options.sarcomere_fields)))
            vi = [vi i];
        end
    end
    % Drop double manipulations
    bi = find( ~strcmp(d.class_2, ""));
    vi = setdiff(vi, bi);
    % Drop excessive heart rates or simulations that stopped early
    bi = find( (d.hr_heart_rate_bpm > options.max_heart_rate) | ...
                (d.time < options.scan_parameters_min_time) );
    vi_contractility = setdiff(vi, bi);

    hold on;
    h_1(1) = scatter(options.x_scaling_factor * d.vent_chamber_radius(vi_contractility), ...
                    options.y_scaling_factor * d.vent_wall_thickness(vi_contractility), ...
                    options.marker_size, ...
                    d.dpdt_max(vi_contractility), ...
                    "filled", ...
                    MarkerFaceAlpha = options.marker_face_alpha, ...
                    MarkerEdgeColor = options.marker_edge_color, ...
                    MarkerEdgeAlpha = options.marker_edge_alpha);

    c = colorbar;
    c.Limits = options.c_limits;
    c.Ticks = options.c_limits([1 end]);
    c.Position = c.Position + [0.11 0.04 0 -0.45];
    c.Label.String = "    dP/dt" + newline + "(mmHg s^{-1})";
    c.Label.HorizontalAlignment = "center";
    c.Label.Rotation = 0;
    c.Label.FontSize = options.legend_font_size;;
    c.Label.Position([1 2]) = c.Label.Position([1 2]) + [0.5 1045];

    % Add in mitochondria
    vi = find( strcmp(d.class_1, "mitochondria") & ...
                (d.total_change_1 < 0) & ...
                (d.time > options.scan_parameters_min_time) & ...
                (d.hr_heart_rate_bpm <= options.max_heart_rate) );

    x = options.x_scaling_factor * d.vent_chamber_radius(vi);
    y = options.y_scaling_factor * d.vent_wall_thickness(vi);

    h_1(2) = scatter(x, y, ...
                options.marker_size_mitochondria, ...
                options.marker_mitochondria, ...
                MarkerFaceColor = options.marker_face_color_mitochondria, ...
                MarkerFaceAlpha = options.marker_face_alpha_mitochondria, ...
                MarkerEdgeColor = options.marker_edge_color_mitochondria);

    p = polyfit(x, y, 1);
    x_fit = linspace(min(x), max(x), 100);
    y_fit = polyval(p, x_fit);
    h_temp = plot(x_fit, y_fit, '-', ...
                Color = options.marker_face_color_mitochondria);
    uistack(h_temp, 'bottom');

    if (options.mybpc_annotation)
        % Highlight c_protein that destabilizes myosin
        vi_c = find(startsWith(d.variable_1, 'c_state_1_1_1'));
        [s, si] = sort(d.total_change_1(vi_c));
        xc(1) = d.vent_chamber_radius(vi_c(si(1)));
        yc(1) = d.vent_wall_thickness(vi_c(si(1)));

        vi_c2 = find(startsWith(d.variable_1, 'c_state_2_1_1'));
        [s, si] = sort(d.total_change_1(vi_c2), 'descend');
        xc(2) = d.vent_chamber_radius(vi_c2(si(1)));
        yc(2) = d.vent_wall_thickness(vi_c2(si(1)));
    end

    % Switch to power
    subplot(subplots(2));

    % Plot reference
    h_temp = scatter(options.x_scaling_factor * d.vent_chamber_radius(vi_contractility), ...
                    options.y_scaling_factor * d.vent_wall_thickness(vi_contractility), ...
                    options.marker_size_contractility, ...
                    options.marker_contractility, ...
                    MarkerFaceColor = options.marker_face_color_contractility, ...
                    MarkerFaceAlpha = options.marker_face_alpha_contractility, ...
                    MarkerEdgeColor = options.marker_edge_color_contractility);
    h_temp.LegendDisplay = "off";

    % Add in blood pressure
    vi = find(startsWith(d.class_1, "baroreflex") & ...
                (d.total_change_1 > 0) & ...
                (d.time > options.scan_parameters_min_time) & ...
                (d.hr_heart_rate_bpm <= options.max_heart_rate) );

    x = options.x_scaling_factor * d.vent_chamber_radius(vi);
    y = options.y_scaling_factor * d.vent_wall_thickness(vi);

    h_2(1) = scatter(x, y, ...
                options.marker_size_hypertension, ...
                options.marker_hypertension, ...
                MarkerFaceColor = options.marker_face_color_hypertension, ...
                MarkerFaceAlpha = options.marker_face_alpha_hypertension, ...
                MarkerEdgeColor = options.marker_edge_color_hypertension);

    p = polyfit(x, y, 1);
    x_fit = linspace(min(x), max(x), 100);
    y_fit = polyval(p, x_fit);
    h_temp = plot(x_fit, y_fit, '-', ...
                Color = options.marker_face_color_hypertension);
    uistack(h_temp, 'bottom');

    % Add in aortic insufficiency
    vi = find(startsWith(d.class_1, "valve") & ...
                strcmp(d.variable_1, "av_valve_leak") & ...
                (d.time > options.scan_parameters_min_time) & ...
                (d.hr_heart_rate_bpm <= options.max_heart_rate) );

    x = options.x_scaling_factor * d.vent_chamber_radius(vi);
    y = options.y_scaling_factor * d.vent_wall_thickness(vi);

    h_2(2) = scatter(x, y, ...
                options.marker_size_av_leak, ...
                options.marker_av_leak, ...
                MarkerFaceColor = options.marker_face_color_av_leak, ...
                MarkerFaceAlpha = options.marker_face_alpha_av_leak, ...
                MarkerEdgeColor = options.marker_edge_color_av_leak);

    p = polyfit(x, y, 1);
    x_fit = linspace(min(x), max(x), 100);
    y_fit = polyval(p, x_fit);
    h_temp = plot(x_fit, y_fit, '-', ...
                Color = options.marker_face_color_av_leak);
    uistack(h_temp, 'bottom');

    % Add in clean mitral insufficiency
    vi = find(startsWith(d.class_3, "valve") & ...
                strcmp(d.variable_3, "mv_valve_leak") & ...
                (d.total_change_1 == 0) & ...
                (d.total_change_2 == 0) & ...
                (d.time > options.scan_parameters_min_time) & ...
                (d.hr_heart_rate_bpm <= options.max_heart_rate) );

    % Plot mv_leak with reduced contractility
    leaks = unique(d.total_change_3(vi));

    for j = 1 : numel(leaks)
        leak_vi = find(startsWith(d.class_3, "valve") & ...
                    strcmp(d.variable_3, "mv_valve_leak") & ...
                    (d.total_change_3 == leaks(j)) & ...
                    (d.time > options.scan_parameters_min_time) & ...
                    (d.hr_heart_rate_bpm <= options.max_heart_rate) );

        base_vi = find( (d.total_change_1(leak_vi) == 0) & ...
                            (d.total_change_2(leak_vi) == 0) );
        base_leak_dpdt = d.dpdt_max(leak_vi(base_vi));

        leak_vi = leak_vi(d.dpdt_max(leak_vi) <= base_leak_dpdt);

        % leak_vi = leak_vi(1 : 2 :end)

        x = options.x_scaling_factor * d.vent_chamber_radius(leak_vi);
        y = options.y_scaling_factor * d.vent_wall_thickness(leak_vi);

        h_2(4) = scatter(x, y, ...
                    options.marker_size_mv_leak_reduced, ...
                    d.dpdt_max(leak_vi), ...
                    "filled", ...
                    options.marker_mv_leak_reduced, ...
                    MarkerFaceAlpha = options.marker_face_alpha_mv_leak_reduced, ...
                    MarkerEdgeColor = options.marker_edge_color_mv_leak_reduced);


        p = polyfit(x, y, 3);
        x_fit = linspace(min(x), max(x), 100);
        y_fit = polyval(p, x_fit);
        h_temp = plot(x_fit, y_fit, '-', ...
            Color = 0.5 * ones(1,3));
        uistack(h_temp, 'bottom');
    end

    x = options.x_scaling_factor * d.vent_chamber_radius(vi);
    y = options.y_scaling_factor * d.vent_wall_thickness(vi);

    h_2(3) = scatter(x, y,  ...
                options.marker_size_mv_leak_clean, ...
                options.marker_mv_leak_clean, ...
                MarkerFaceColor = options.marker_face_color_mv_leak_clean, ...
                MarkerEdgeColor = options.marker_edge_color_mv_leak_clean);

    p = polyfit(x, y, 1);
    x_fit = linspace(min(x), max(x)+0.5, 100);
    y_fit = polyval(p, x_fit);
    h_temp = plot(x_fit, y_fit, '-', ...
        Color = options.marker_face_color_mv_leak_clean);
    uistack(h_temp, 'bottom');


    c2 = colorbar;
    c2.Limits = options.c_limits_2;
    c2.Ticks = options.c_limits_2([1 end]);
    c2.FontSize = options.legend_font_size;
    c2.Position = c2.Position + [0.11 0.04 0 -0.45];
    c2.Label.String = "    dP/dt" + newline + "(mmHg s^{-1})";
    c2.Label.HorizontalAlignment = "center";
    c2.Label.Rotation = 0;
    c2.Label.FontSize = options.legend_font_size;;
    c2.Label.Position([1 2]) = c2.Label.Position([1 2]) + [0.5 700];

    % Set the legends
    subplot(subplots(1));
    legend_uktcvr(h_1(1:2), options.lab_1, ...
            legend_position = options.legend_position, ...
            legend_vertical_scale = options.legend_vertical_scale, ...
            IconColumnWidth = options.legend_icon_column_width, ...
            FontSize = options.legend_font_size);

    subplot(subplots(2));
    legend_uktcvr(h_2, options.lab_2, ...
            legend_position = options.legend_position_2, ...
            legend_vertical_scale = options.legend_vertical_scale_2, ...
            IconColumnWidth = options.legend_icon_column_width, ...
            FontSize = options.legend_font_size)


    for i = 1 : 2
        ax = improve_axes( ...
            'axis_handle', subplots(i), ...
            'x_axis_label', 'Diastolic chamber radius (mm)', ...
            'x_label_offset', options.x_label_offset, ...
            'x_ticks', options.x_ticks, ...
            'x_tick_decimal_places', 0, ...
            'x_tick_length', 0.025, ...
            'y_axis_offset', 0.1, ...
            'y_axis_label', {'Diastolic', 'wall','thickness', '(mm)'}, ...
            'y_label_offset', options.y_label_offset, ...
            'y_ticks', options.y_ticks, ...
            'y_tick_decimal_places', 0, ...
            'tick_font_size', options.tick_font_size, ...
            'title', options.title_strings{i}, ...
            'title_y_offset', options.title_y_offset);
    end
    
    

    % Add in annotations

    % Base model
    subplot(subplots(1));

    arrow_start = [options.base_annotation_x options.base_annotation_y];
    ha = arrow(arrow_start, ...
                [options.x_scaling_factor * base_radius ...
                    options.y_scaling_factor * base_thickness], ...
                'length', 10, ...
                'base', 50);

    ha.FaceColor = options.base_annotation_color;
    ha.EdgeColor = options.base_annotation_color;
    ha.LegendDisplay = 'off';

    t = text(arrow_start(1) + options.base_annotation_text_x_offset, ...
            arrow_start(2) + options.base_annotation_text_y_offset, ...
            options.base_annotation_label, ...
            FontSize = options.annotation_font_size, ...
            HorizontalAlignment = 'center', ...
            VerticalAlignment = 'bottom', ...
            Color = options.base_annotation_color);

    if (options.mybpc_annotation)

        subplot(subplots(1))

        for ci = 1 : numel(xc)
            arrow_start = [options.mybpc_annotation_x options.mybpc_annotation_y];
            ha = arrow(arrow_start, ...
                [options.x_scaling_factor * xc(ci) options.y_scaling_factor * yc(ci)], ...
                'length', 10, ...
                'base', 50);

            ha.FaceColor = options.mybpc_annotation_color;
            ha.EdgeColor = options.mybpc_annotation_color;
            ha.LegendDisplay = 'off';
        end

        t = text(arrow_start(1) + options.mybpc_annotation_text_x_offset, ...
            arrow_start(2) + options.mybpc_annotation_text_y_offset, ...
            options.mybpc_annotation_label, ...
            FontSize = options.annotation_font_size, ...
            HorizontalAlignment = 'center', ...
            VerticalAlignment = 'bottom', ...
            Color = options.mybpc_annotation_color);
    end

     % Output figure
    if (~(options.output_image_file == ""))
        for out_i = 1 : numel(options.output_image_types)
            output_file_string = ...
                    sprintf("%s.%s", ...
                        options.output_image_file, ...
                        options.output_image_types(out_i));
            
            disp(sprintf('Writing image to: %s', output_file_string));
    
            exportgraphics(fig_options.figure_handle, ...
                output_file_string, ...
                Resolution = 1200, ...
                ContentType = "auto");
        end
    end