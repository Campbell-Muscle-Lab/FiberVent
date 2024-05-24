# FiberVent

````mermaid

classDiagram

    CMyoVent <|-- cmv_system
    cmv_system <|-- cmv_model
    cmv_system <|-- cmv_protocol
    cmv_system <|-- cmv_results
    cmv_system <|-- circulation
    circulation <|-- hemi_vent
    circulation <|-- VAD
    circulation <|-- baroreflex
    circulation <|-- Growth
    hemi_vent <|-- muscle
    hemi_vent <|-- valve
    muscle <|-- half_sarcomere
    half_sarcomere <|-- membranes
    half_sarcomere <|-- heart_rate
    half_sarcomere <|-- myofilaments
    half_sarcomere <|-- mitochondria
    myofilaments <|--kinetic_scheme
    kinetic_scheme <|-- m_state
    m_state <|-- transition
    baroreflex <|-- reflex_control
    Growth <|-- Growth_control
    cmv_protocol <|-- activation
    cmv_protocol <|-- perturbation
    
    %% cmv_system: + *cmv_model p_cmv_model
    %% cmv_system: + run_simulation()

    %% circulation: +circ_blood_volume
    %% circulation: +circ_no_of_compartments
    %% circulation: +circ_resistance[]
    %% circulation: +circ_compliance[]
    %% circulation: +circ_slack_volume[]
    %% circulation: +circ_inertance[]
    %% circulation: +circ_pressure[]
    %% circulation: +circ_volume[]
    %% circulation: +circ_flow[]
    %% circulation: +circ_total_slack_volume

    %% cmv_results: +no_of_time_points
    %% cmv_results: +no_of_beats
    %% cmv_results: +last_beat_t_index
    %% cmv_results: +field_indices_[many]

    %% hemi_vent: +vent_wall_density
    %% hemi_vent: +vent_wall_volume
    %% hemi_vent: +vent_chamber_volume
    %% hemi_vent: +vent_chamber_pressure
    %% hemi_vent: +vent_wall_thickness
    %% hemi_vent: +vent_n_hs
    %% hemi_vent: +vent_circumference
    %% hemi_vent: +vent_stroke_work_J
    %% hemi_vent: +vent_energy_used_J
    %% hemi_vent: +vent_efficiency
    %% hemi_vent: +vent_ATP_used_per_s

    %% half_sarcomere: + hs_prop_fibrosis
    %% half_sarcomere: + hs_prop_myofilaments
    %% half_sarcomere: + hs_length
    %% half_sarcomere: + hs_stress
    %% half_sarcomere: + hs_ATP_concentration
    %% half_sarcomere: + hs_delta_G_ATP

    %% valve: + valve_pos
    %% valve: + valve_vel
    %% valve: + valve_name
    %% valve: + valve_mass
    %% valve: + valve_eta
    %% valve: + valve_k

    %% mitochondria: + volume
    %% mitochondria: + mito_ATP_generation_rate
    %% mitochondria: + mito_ATP_generated_M_per_s



````
