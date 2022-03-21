clear all;
close all;

MATERIALS = ["steel_vase", "kitchen_sponge", "flour_sack", "car_sponge", "black_foam", "acrylic"];
COLOURS = ["red", "green", "blue", "cyan", "magenta", "yellow"];
MARKER_TYPES = ['o', '+', '*', 's', 'd', '^'];
PCA_COLOURS = ["red", "green", "blue"];

STEEL_VASE = dir("steel_vase*.mat");
KITCHEN_SPONGE = dir("kitchen_sponge*.mat");
FLOUR_SACK = dir("flour_sack*.mat");
CAR_SPONGE = dir("car_sponge*.mat");
BLACK_FOAM = dir("black_foam*.mat");
ACRYLIC = dir("acrylic*.mat");

SEPARATE_MATERIALS = [STEEL_VASE, KITCHEN_SPONGE, FLOUR_SACK, CAR_SPONGE, BLACK_FOAM, ACRYLIC];

% "F0Electrodes","F1Electrodes", - Electrode Impedance
% "F0pac","F1pac",                   - High Frequency Fluid Vibrations
% "F0pdc","F1pdc",                   - Low Frequency Fluid Pressure
% "F0tac","F1tac",                   - Core Temperature Change
% "F0tdc","F1tdc",                   - Core Temperature
% "JEff",                        - Robot arm joint effort (load)
% "JPos"  – Robot arm joint positions
% "JVel"                         - Robot arm joint velocity

% The Pac variable is 22-dimensional, but should be 1-dimensional. Please only use the second row when sampling. Thanks to Ezgi for spotting this.  

%% Section A: Data Preparation - [10 marks]

% 1. Use the plot command to view the time series sensor data for the variables Pressure, Vibration
% and Temperature (PVT) and the Electrodes. Do this for several objects and trials and then
% choose a single time step that looks like it will allow differentiation between the data for
% different objects. Explain why you chose that value. Include an example of your data
% visualisation for one or two object trials in your report.

% Plotting all trials for all materials for each variable

plotallvar("Pressure", 1, MATERIALS);
plotallvar("Vibrations", 2, MATERIALS);
plotallvar("Temperature", 3, MATERIALS);

TRIAL_NUM = 5;

% Overlaying P, V and T data for different materials and trial TRIAL_NUM

fig = figure(4);
trial_subplot(fig, 1, TRIAL_NUM, SEPARATE_MATERIALS, "Pressure", 1, 1000);
trial_subplot(fig, 2, TRIAL_NUM, SEPARATE_MATERIALS, "Vibrations", 1, 1000);
trial_subplot(fig, 3, TRIAL_NUM, SEPARATE_MATERIALS, "Temperature", 1, 1000);

% Electrode data for different materials and trial TRIAL_NUM

plot_electrodes(5, TRIAL_NUM, SEPARATE_MATERIALS, MATERIALS);

% Overlaying P, V and T data for different materials and trial TRIAL_NUM at
% time step TIME_STEP

TIME_STEP = 400;

fig = figure(6);
trial_subplot(fig, 1, TRIAL_NUM, SEPARATE_MATERIALS, "Pressure", TIME_STEP, TIME_STEP+1);
trial_subplot(fig, 2, TRIAL_NUM, SEPARATE_MATERIALS, "Vibrations", TIME_STEP, TIME_STEP+1);
trial_subplot(fig, 3, TRIAL_NUM, SEPARATE_MATERIALS, "Temperature", TIME_STEP, TIME_STEP+1);

%%

% 2. For one finger (F0 or F1), sample the Pressure, Vibration, Temperature time series data into
% scaler values measured at the time instance (of your selected time step) for each object / trial.
% Save the data structures together as a .mat file called F0_PVT.mat or F1_PVT.mat. Repeat for
% the Electrodes data, saving that as another .mat file. Note that all subsequent actions in this
% coursework will be on the data sets you just created (and therefore only on one of the robot’s
% fingers).

F1_PVT = zeros(6, 10, 3);  % # materials x # trials x # variables
F1_ELECTRODES = zeros(6, 10, 19); % # materials x # trials x # variables

for i = 1:6
    material = SEPARATE_MATERIALS(:, i);
    for j = 1:10
        data = load(material(j).name);
        F1_PVT(i,j,1) = data.F1pdc(1,TIME_STEP);  % Pressure
        F1_PVT(i,j,2) = data.F1tdc(1,TIME_STEP);  % Temperature
        F1_PVT(i,j,3) = data.F1pac(2,TIME_STEP);  % Vibrations
    end
end

save("F1_PVT.mat", "F1_PVT");

for i = 1:6
    material = SEPARATE_MATERIALS(:, i);
    for j = 1:10
        data = load(material(j).name);
        F1_ELECTRODES(i,j,:) = data.F1Electrodes(:,TIME_STEP);  % Electrodes
    end
end

save("F1_ELECTRODES.mat", "F1_ELECTRODES");

% 3. Create a 3D scatter plot of the complete contents of the PVT mat file, with the axis as Pressure,
% Vibration and Temperature, with different colours used for different objects. Use the same
% colours for the objects throughout this work.

fig = figure(7);
ax = subplot(1, 1, 1, "Parent", fig);
for i = 1:10
    for j = 1:6
        scatter3(ax, F1_PVT(j, i, 1), F1_PVT(j, i, 2), F1_PVT(j, i, 3), "o", COLOURS(j), "filled");
        hold(ax, "on");
    end
end
xlabel(ax, "Pressure");
ylabel(ax, "Vibration");
zlabel(ax, "Temperature");
title(ax, "PVT Scatter Plot");
legend(ax, "Steel Vase", "Kitchen Sponge", "Flour Sack", "Car Sponge", "Black Foam", "Acrylic");
hold(ax, "off");

%% Section B: Principal Component Analysis – [25 marks]

% 1. Using PCA (Principal Component Analysis) determine the principal components of the PVT
% data.
% a. Report covariance matrix, eigenvalues, and eigenvectors for the data.
% b. Replot the Standardised data with the Principal components displayed.
% c. Reduce the data to 2-dimensions and replot.
% d. Show how the data is distributed across all principal components by plotting as
% separate 1D number lines.
% e. Comment on your findings.

% Reshape data by removing class dimension
pvt_3d_points = reshape(F1_PVT, [60, 3]);

% Standardise data by subtracting mean and dividing by standard deviation
pvt_mean = mean(pvt_3d_points, 1);
pvt_std = std(pvt_3d_points, 1);
pvt_3d_points_normalised = (pvt_3d_points - pvt_mean) ./ pvt_std;

% Calculate covariance matrix
C = cov(pvt_3d_points_normalised);
disp(C);

% Calculate eigen vectors/values of covariance matrix
[V,D] = eig(C);
disp(V);
disp(D);

% Sort eigenvectors of covariance matrix in order of decreasing eigenvalue
% magnitude 
[~,idx] = sort(diag(D), 'descend');
pca_coeff = V(:, idx);

% Prepare PCA data from plotting vectors as lines from origin
lines = zeros(2,3,3);
lines(2, :, :) = pca_coeff;

% Reshape PVT data so that classes can be colour coded
pvt_3d_points_normalised_reshaped = reshape(pvt_3d_points_normalised, size(F1_PVT));

% Create scatter plot
fig = figure(8);
ax = subplot(1, 1, 1, "Parent", fig);
for i = 1:6
    scatter3(ax, pvt_3d_points_normalised_reshaped(i, :, 1), pvt_3d_points_normalised_reshaped(i, :, 2), pvt_3d_points_normalised_reshaped(i, :, 3), "o", COLOURS(i), "filled");
    hold(ax, "on");
end
plot3(ax, lines(:,:, 1), lines(:,:, 2), lines(:,:, 3), "linewidth", 2);
xlabel(ax, "Pressure");
ylabel(ax, "Vibration");
zlabel(ax, "Temperature");
title(ax, "Standardised PVT Scatter Plot with principle components");
legend(ax, "Steel Vase", "Kitchen Sponge", "Flour Sack", "Car Sponge", "Black Foam", "Acrylic", "PC1", "PC2", "PC3");
hold(ax, "off");

% Project data onto principle components
projected = pvt_3d_points_normalised * pca_coeff;

% Show scatter plot of data projected onto top 2 PC's
projected_2d = projected(:, 1:2);
projected_2d_reshaped = reshape(projected_2d, [6, 10, 2]);
fig = figure(9);
ax = subplot(1, 1, 1, "Parent", fig);
for i = 1:6
    scatter(ax, projected_2d_reshaped(i, :, 1), projected_2d_reshaped(i, :, 2), "o", COLOURS(i), "filled");
    hold(ax, "on");
end
xlabel(ax, "PC1");
ylabel(ax, "PC2");
title(ax, "Scatter Plot of PVT data projected down with 2 principle components");
legend(ax, "Steel Vase", "Kitchen Sponge", "Flour Sack", "Car Sponge", "Black Foam");
hold(ax, "off");

% Show 3 1D plot of data projected onto each PC
projected_reshaped = reshape(projected, [6, 10, 3]);
fig = figure(10);
plot_pc(fig, 1, projected_reshaped(:,:,1), COLOURS);
plot_pc(fig, 2, projected_reshaped(:,:,2), COLOURS);
ax = plot_pc(fig, 3, projected_reshaped(:,:,3), COLOURS);
legend(ax, "Steel Vase", "Kitchen Sponge", "Flour Sack", "Car Sponge", "Black Foam", "Acrylic");

%%

% 2. There are 19 electrodes per sensor, so relationship between the electrodes for different
% objects cannot be easily visualised as in the last questions.
% a. Use PCA to determine the principal components of the electrode data. Report on the
% variances of each principal components using a Scree plot.
% b. Visualize the electrode data using the three principal components with largest
% variance.
% c. Comment on your findings. 

% Reshape data by removing class dimension
electrode_3d_points = reshape(F1_ELECTRODES, [60, 19]);

% Standardise data by subtracting mean and dividing by standard deviation
electrode_mean = mean(electrode_3d_points, 1);
electrode_std = std(electrode_3d_points, 1);
electrode_3d_points_normalised = (electrode_3d_points - electrode_mean) ./ electrode_std;

% Calculate covariance matrix
C = cov(electrode_3d_points_normalised);

% Calculate eigen vectors/values of covariance matrix
[V,D] = eig(C);

% Sort eigenvectors of covariance matrix in order of decreasing eigenvalue
% magnitude 
[eigs,idx] = sort(diag(D), 'descend');
electrode_pca_coeff = V(:, idx);

% Show variances of each principal component using a Scree plot
fig = figure(11);
ax = subplot(1, 1, 1, "Parent", fig);
plot(ax, 1:size(eigs), eigs, '-o');
title(ax, "Scree plot");
xlabel(ax, "Principle component");
ylabel(ax, "Variance");
xticks(ax, 1:size(eigs));

% Visualize the electrode data using the three principal components with largest variance.
electrodes_projected = reshape(electrode_3d_points_normalised * electrode_pca_coeff, size(F1_ELECTRODES));
fig = figure(12);
ax = subplot(1, 1, 1, "Parent", fig);
for i = 1:6
    scatter3(ax, electrodes_projected(i, :, 1), electrodes_projected(i, :, 2), electrodes_projected(i, :, 3), "o", COLOURS(i), "filled");
    hold(ax, "on");
end
xlabel(ax, "PC 1");
ylabel(ax, "PC 2");
zlabel(ax, "PC 3");
title(ax, "Electrode data projected onto top 3 principle components");
legend(ax, "Steel Vase", "Kitchen Sponge", "Flour Sack", "Car Sponge", "Black Foam", "Acrylic");
hold(ax, "off");


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

% K means on normalised data, with k = 6 clusters, different distances
idx1 = clustering(pvt_3d_points_normalised, COLOURS, MARKER_TYPES, 'sqeuclidean');
disp(clustering_accuracy(idx1));
idx2 = clustering(pvt_3d_points_normalised, COLOURS, MARKER_TYPES, 'cityblock');
disp(clustering_accuracy(idx2));
idx3  = clustering(pvt_3d_points_normalised, COLOURS, MARKER_TYPES, 'cosine');
disp(clustering_accuracy(idx3));
idx4 = clustering(pvt_3d_points_normalised, COLOURS, MARKER_TYPES, 'correlation');
disp(clustering_accuracy(idx4));

% 2. Now apply bagging (bootstrap aggregation for an ensemble of decision trees) to the electrode
% data that was processed with PCA in section B.2.b. Use a 60 / 40 split for Training / Test data.
% a. Specify the number of bags / trees you used. Why did you choose this number?
% b. Visualise two of your generated decision trees.
% c. Run the trained model with the test data. Display a confusion matrix (where the
% object type is the class) and comment on the overall accuracy.
% d. Discuss the following: How can misclassifications in your results be explained given
% the object properties? Do you think the PCA step was helpful?

% Split into train and test set
[m,n] = size(projected) ;
split = 0.60;
idx = randperm(m)  ;
train = projected(idx(1:round(split*m)),:); 
test = projected(idx(round(split*m)+1:end),:);

% Make labels
labels = ones(6, 10, 1);
for i = 1:6
    labels(i,:,:) = i;
end
labels = reshape(labels, 60, 1);
train_labels = labels(idx(1:round(split*m)),:);
test_labels = labels(idx(round(split*m)+1:end),:);

% Train a bagging moddel
NUM_TREES = 12;
B = TreeBagger(NUM_TREES,train,train_labels, 'OOBPrediction', 'on');
view(B.Trees{1},'Mode','graph');
view(B.Trees{2},'Mode','graph');

% Plot OOB error to find optimal number of trees (works when NUM_TREES high)
% Shows a local minima at 12
figure();
oobErrorBaggedEnsemble = oobError(B);
plot(oobErrorBaggedEnsemble, "linewidth", 2);
xticks(1:12);
xlabel('Number of grown trees');
ylabel('Out-of-bag classification error');

% Test model and display confusion matrix
preds = cellfun(@str2num,predict(B, test));
C = confusionmat(test_labels,preds);
confusionchart(C);

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
% recognition. What are the pros and cons of this other approach


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
                    xlabel(ax, "Time (ms)");
                    ylabel(ax, "Pressure");
                end
                if varname == "Temperature"
                    [~, n] = size(data.F1tdc);
                    plot(ax, 1:n, data.F1tdc(1,:));
                    xlabel(ax, "Time (ms)");
                    ylabel(ax, "Temperature");
                end
                if varname == "Vibrations"
                    [~, n] = size(data.F1pac);
                    plot(ax, 1:n, data.F1pac(2,:));
                    xlabel(ax, "Time (ms)");
                    ylabel(ax, "Vibration");
                end
                title(ax, varname+" "+ strrep(material, "_", " "));
                hold(ax, "on");
            end
        end
    end
    legend(ax, "Trial 1", "Trial 2", "Trial 3", "Trial 4", "Trial 5", "Trial 6", "Trial 7", "Trial 8", "Trial 9", "Trial 10");
end

function ax = trial_subplot(fig, sub_plot_num, trial_num, seperate_materials, varname, t_start, t_end)
    ax = subplot(3, 1, sub_plot_num, "Parent", fig);
    hold(ax, "on");
    for i = 1:6
        file = seperate_materials(trial_num, i).name;
        data = load(file);
        if varname == "Pressure"
            plot(ax, t_start:t_end, data.F1pdc(1,t_start:t_end));
            xlabel(ax, "Time (ms)");
            ylabel(ax, "Pressure");
        end
        if varname == "Temperature"
            plot(ax, t_start:t_end, data.F1tdc(1,t_start:t_end));
            xlabel(ax, "Time (ms)");
            ylabel(ax, "Temperature");
        end
        if varname == "Vibrations"
            plot(ax, t_start:t_end, data.F1pac(2,t_start:t_end));
            xlabel(ax, "Time (ms)");
            ylabel(ax, "Vibration");
        end
    end
    xticks(ax, [t_start, t_end]);
    title(ax, varname);
    legend(ax, "Steel Vase", "Kitchen Sponge", "Flour Sack", "Car Sponge", "Black Foam", "Acrylic");
    hold(ax, "off");
end

function fig = plot_electrodes(fignum, trial_num, seperate_materials, material_names)
    fig = figure(fignum);
    for i = 1:6
        file = seperate_materials(trial_num, i).name;
        data = load(file);
        ax = subplot(2,3, i, "Parent", fig);
        elec_data = data.F1Electrodes;
        [n_electrodes, n_timesteps] = size(elec_data);
        for j = 1:n_electrodes
            plot(ax, 1:n_timesteps, elec_data(j,:));
            hold(ax, "on");
        end
        hold(ax, "off");
        title(ax, strrep(material_names(1, i), "_", " "));
        xlabel(ax, "Time (ms)");
        ylabel(ax, "Impedance");
    end
    legend("Electrode 1","Electrode 2","Electrode 3","Electrode 4","Electrode 5","Electrode 6","Electrode 7","Electrode 8","Electrode 9","Electrode 10","Electrode 11","Electrode 12","Electrode 13","Electrode 14","Electrode 15","Electrode 16","Electrode 17","Electrode 18","Electrode 19");
end

function ax = plot_pc(fig, pc_num, pc_data, colours)
    ax = subplot(3,1, pc_num, "Parent", fig);
    for i = 1:6
        scatter(ax, pc_data(i,:), zeros(size(pc_data(i,:))), "o", colours(i), "filled");
        hold(ax, "on");
    end
    title(ax, "PC "+num2str(pc_num));
    xlim([-5, 5]);
end

function idx = clustering(pvt_3d_points_normalised, colours, marker_types, distance_metric)
    idx = kmeans(pvt_3d_points_normalised, 6, "Distance", distance_metric); 

    fig = figure();
    ax = subplot(1, 1, 1, "Parent", fig);

    % Plot ground truth labels
    reshaped = reshape(pvt_3d_points_normalised, [6, 10, 3]);
    for j = 1:6
        scatter3(ax, reshaped(j, :, 1), reshaped(j, :, 2), reshaped(j, :, 3), 80, "o", colours(j), "filled");
        hold(ax, "on");
    end

    % Plot clusters
    for j = 1:6
        data = pvt_3d_points_normalised(idx == j, :);
        scatter3(ax, data(:, 1), data(:, 2), data(:, 3), 40, marker_types(j), "black", "linewidth", 1.5);
    end
    legend(ax, "Steel Vase", "Kitchen Sponge", "Flour Sack", "Car Sponge", "Black Foam", "Acrylic", "1", "2", "3", "4", "5", "6");
    xlabel(ax, "Pressure");
    ylabel(ax, "Vibration");
    zlabel(ax, "Temperature");
    title(ax, "Clustered data with k=6 clusters, "+distance_metric+" distance metric");
    hold(ax, "off");
end

function acc = clustering_accuracy(idx)
    labels1 = [1,2,3,4,5,6]';
    labels2 = [2,3,4,5,6,1]';
    labels3 = [3,4,5,6,1,2]';
    labels4 = [4,5,6,1,2,3]';
    labels5 = [5,6,1,2,3,4]';
    labels6 = [2,3,4,5,6,1]';
    l1 = repmat(labels1,10,1);
    l2 = repmat(labels2,10,1);
    l3 = repmat(labels3,10,1);
    l4 = repmat(labels4,10,1);
    l5 = repmat(labels5,10,1);
    l6 = repmat(labels6,10,1);
    acc1 = sum(idx == l1) / 60;
    acc2 = sum(idx == l2) / 60;
    acc3 = sum(idx == l3) / 60;
    acc4 = sum(idx == l4) / 60;
    acc5 = sum(idx == l5) / 60;
    acc6 = sum(idx == l6) / 60;
    acc = max([acc1, acc2, acc3, acc4, acc5, acc6]);
end