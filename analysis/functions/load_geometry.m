function [bodies, num_bodies, rotation_face_index, x_label] = load_geometry(geometry_type, energy_accommodation, surface_temperature__K, show_geometry)
%LOAD_GEOMETRY loads a geometry based on the specified type
%   [bodies, num_bodies, rotation_face_index, x_label] = LOAD_GEOMETRY(geometry_type, energy_accommodation, surface_temperature__K, show_geometry)
%   loads the geometry specified by geometry_type and returns the bodies,
%   number of bodies, rotation face index, and x-axis label.
%   If show_geometry is true, it will display the geometry.

    if nargin < 4
        show_geometry = true;
    end

    if strcmp(geometry_type, 'plate')
        bodies = parametrized_flat_plate(1, 1, [0,0,0],true,energy_accommodation,surface_temperature__K);
        if show_geometry
            showBodies(bodies, [0], 0.75, 0.25);
        end
        num_bodies = 1;
        rotation_face_index = 1;
        x_label = "angle of attack [°]";
    elseif strcmp(geometry_type, 'shuttlecock')
        bodies = load_from_gmsh(energy_accommodation,surface_temperature__K);
        if show_geometry
            showBodies(bodies, [0,pi/4,pi/4,pi/4,pi/4], 0.75, 0.25);
        end
        num_bodies = 5;
        rotation_face_index = [2,3,4,5];
        x_label = "control surface angle [°]";
    elseif strcmp(geometry_type, 'shuttlecock_wing')
        bodies = load_shuttlecock_wing(energy_accommodation,surface_temperature__K);
        if show_geometry
            showBodies(bodies, [0/4], 0.75, 0.25);
        end
        num_bodies = 1;
        rotation_face_index = 1;
        x_label = "angle of attack [°]";
    elseif strcmp(geometry_type, 'shuttlecock_wing_new')
        bodies_all = load_from_gmsh(energy_accommodation,surface_temperature__K);
        bodies = cell(1,1);
        bodies{1} = bodies_all{3};
        if show_geometry
            showBodies(bodies, [pi/4], 0.75, 0.25);
        end
        num_bodies = 1;
        rotation_face_index = 1;
        x_label = "angle of attack [°]";
    elseif strcmp(geometry_type,'box')
        bodies = load_box(energy_accommodation,surface_temperature__K);
        if show_geometry
            showBodies(bodies,[0]);
        end
        num_bodies = 1;
        x_label = "angle of attack [°]";
        rotation_face_index = 1;
    else
        error("Invalid geometry_type. Use 'plate', 'shuttlecock', 'shuttlecock_wing', 'shuttlecock_wing_new' or 'box'.");
    end
end
