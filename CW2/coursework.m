MATERIALS = ["steel_vase", "kitchen_sponge", "flour_sack", "car_sponge", "black_foam", "acrylic"];
STEEL_VASE = dir("steel_vase*.mat");
KITCHEN_SPONGE = dir("kitchen_sponge*.mat");
FLOUR_SACK = dir("flour_sack*.mat");
CAR_SPONGE = dir("car_sponge*.mat");
BLACK_FOAM = dir("black_foam*.mat");
ACRYLIC = dir("acrylic*.mat");

SEPARATE_MATERIALS = [STEEL_VASE, KITCHEN_SPONGE, FLOUR_SACK, CAR_SPONGE, BLACK_FOAM, ACRYLIC];

% 'F0Electrodes','F1Electrodes', - Electrode Impedance
% 'F0pac','F1pac',                   - High Frequency Fluid Vibrations
% 'F0pdc','F1pdc',                   - Low Frequency Fluid Pressure
% 'F0tac','F1tac',                   - Core Temperature Change
% 'F0tdc','F1tdc',                   - Core Temperature
% 'JEff',                        - Robot arm joint effort (load)
% 'JPos'  – Robot arm joint positions
% 'JVel'                         - Robot arm joint velocity

% The Pac variable is 22-dimensional, but should be 1-dimensional. Please only use the second row when sampling. Thanks to Ezgi for spotting this.  

%% Section A: Data Preparation - [10 marks]

% 1. Use the plot command to view the time series sensor data for the variables Pressure, Vibration
% and Temperature (PVT) and the Electrodes. Do this for several objects and trials and then
% choose a single time step that looks like it will allow differentiation between the data for
% different objects. Explain why you chose that value. Include an example of your data
% visualisation for one or two object trials in your report.

plotallvar("Pressure", 1, MATERIALS);
plotallvar("Vibrations", 2, MATERIALS);
plotallvar("Temperature", 3, MATERIALS);

fig = figure(4);
trial_subplot(fig, 1, 1, SEPARATE_MATERIALS, "Pressure", 1, 1000);
trial_subplot(fig, 2, 1, SEPARATE_MATERIALS, "Vibrations", 1, 1000);
trial_subplot(fig, 3, 1, SEPARATE_MATERIALS, "Temperature", 1, 1000);

% 2. For one finger (F0 or F1), sample the Pressure, Vibration, Temperature time series data into
% scaler values measured at the time instance (of your selected time step) for each object / trial.
% Save the data structures together as a .mat file called F0_PVT.mat or F1_PVT.mat. Repeat for
% the Electrodes data, saving that as another .mat file. Note that all subsequent actions in this
% coursework will be on the data sets you just created (and therefore only on one of the robot’s
% fingers).

% 3. Create a 3D scatter plot of the complete contents of the PVT mat file, with the axis as Pressure,
% Vibration and Temperature, with different colours used for different objects. Use the same
% colours for the objects throughout this work.

%% Section B: Principal Component Analysis – [25 marks]

% 1. Using PCA (Principal Component Analysis) determine the principal components of the PVT
% data.
% a. Report covariance matrix, eigenvalues, and eigenvectors for the data.
% b. Replot the Standardised data with the Principal components displayed.
% c. Reduce the data to 2-dimensions and replot.
% d. Show how the data is distributed across all principal components by plotting as
% separate 1D number lines.
% e. Comment on your findings.

% 2. There are 19 electrodes per sensor, so relationship between the electrodes for different
% objects cannot be easily visualised as in the last questions.
% a. Use PCA to determine the principal components of the electrode data. Report on the
% variances of each principal components using a Scree plot.
% b. Visualize the electrode data using the three principal components with largest
% variance.
% c. Comment on your findings. 

%% Section C: Linear Discriminant Analysis (LDA) - [20 marks]

% 1. We want to see if we can discriminate two deformable and porous objects by touch: the black
% foam and the car sponge.
% a. Use LDA to split the training data in terms of Pressure vs. Vibration, Pressure vs.
% Temperature and Temperature vs. Vibration. Plot the results, including a line showing
% the generated LDA function.
% b. Now apply LDA to the three-dimensional PVT data.
% c. Comment on the different outcomes. Consider the physical properties of the objects
% in your answer and how these may have affected the sensor readings.
% d. Repeat the LDA analysis with your own choice of two objects. Explain why you have
% selected those objects for analysis. In other words – what were you trying to test and
% what did you determine?

%% Section D: Clustering & Classification - [30 marks]

% 1. Apply your choice of a clustering algorithm (that we covered in class) to the PVT data
% a. Visualise the output.
% b. Comment on the outcome. Do the clusters correspond to real-life object similarities?
% c. Change the distance metric, repeat the clustering and comment on the change in the
% outcome.

% 2. Now apply bagging (bootstrap aggregation for an ensemble of decision trees) to the electrode
% data that was processed with PCA in section B.2.b. Use a 60 / 40 split for Training / Test data.
% a. Specify the number of bags / trees you used. Why did you choose this number?
% b. Visualise two of your generated decision trees.
% c. Run the trained model with the test data. Display a confusion matrix (where the
% object type is the class) and comment on the overall accuracy.
% d. Discuss the following: How can misclassifications in your results be explained given
% the object properties? Do you think the PCA step was helpful?

%% Section E: Conclusion – [15 marks]

% 1. Summarise your work.
% a. How have the pattern recognition techniques helped us analyse the data?
% b. Would you say it is possible to distinguish the objects only using touch?
% c. If you wanted to repeat the experiment with a cheaper tactile sensor (one with fewer
% sensing modalities and/or electrode channels than the BioTac), what object
% properties do you think would be most important for that sensor to measure? Justify
% your answer based on your findings.
% d. Our analysis is based on a single time step access all the available data sensor data.
% Discuss an alternative method we could use to prepare the data for pattern
% recognition. What are the pros and cons of this other approach?

%% Functions

function fig = plotallvar(varname, fignum, materials)
    [~, n_materials] = size(materials);
    fig = figure(fignum);
    all_files = dir("*.mat");
    for i = 1:size(all_files)
        file = all_files(i);
        for j = 1:n_materials
            material = materials(j);
            if contains(file.name, material)
                data = load(file.name);
                ax = subplot(2,3, j, "Parent", fig);
                if varname == "Pressure"
                    [~, n] = size(data.F1pdc);
                    plot(ax, 1:n, data.F1pdc(1,:));
                end
                if varname == "Temperature"
                    [~, n] = size(data.F1tdc);
                    plot(ax, 1:n, data.F1tdc(1,:));
                end
                if varname == "Vibrations"
                    [~, n] = size(data.F1pac);
                    plot(ax, 1:n, data.F1pac(2,:));
                end
                title(ax, varname+" "+ strrep(material, "_", " "));
                hold(ax, "on");
            end
        end
    end
end

function ax = trial_subplot(fig, sub_plot_num, trial_num, seperate_materials, varname, t_start, t_end)
    ax = subplot(3, 1, sub_plot_num, "Parent", fig);
    hold(ax, "on");
    for i = 1:6
        material = seperate_materials(i);
        file = material(1,trial_num).name;
        data = load(file);
        if varname == "Pressure"
            [~, n] = size(data.F1pdc);
            plot(ax, 1:n, data.F1pdc(1,t_start:t_end));
        end
        if varname == "Temperature"
            [~, n] = size(data.F1tdc);
            plot(ax, 1:n, data.F1tdc(1,t_start:t_end));
        end
        if varname == "Vibrations"
            [~, n] = size(data.F1pac);
            plot(ax, 1:n, data.F1pac(2,t_start:t_end));
        end
    end
    title(ax, varname);
    legend(ax, "Steel Vase", "Kitchen Sponge", "Flour Sack", "Car Sponge", "Black Foam", "Acrylic");
    hold(ax, "off");
end
