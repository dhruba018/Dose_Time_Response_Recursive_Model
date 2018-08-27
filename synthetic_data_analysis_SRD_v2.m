%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Synthetic data example demonstrating the application of Dose-time Recurisve Hybrid Model. 
% The data given in 'syntheic_recursive.mat' has 7 subjects, 7 doses, 8 time points and 14 proteins.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;    clear;      
% close all

% Read data as tables...
PATH = [pwd, '\'];
FILENAME = 'syntheic_recursive.mat';
syn_data = load([PATH, FILENAME]);                 % Read as a structure
dose_time_px = syn_data.pxdata;
dose_time_resp = syn_data.ydata;
time_slope_px = syn_data.beta;

% Parameters...
subjects = unique(dose_time_px.Indiv);                                     m = numel(subjects);
time_pts = unique(dose_time_px.Time_point);                            T = numel(time_pts);
tp_diff = join([time_pts(2:end), time_pts(1:end-1)], '-');
doses = unique(dose_time_px.Dose);                                          D = numel(doses);
proteins = dose_time_px.Properties.VariableNames(4:end)';       P = numel(proteins);


%%
% Prediction for a specific time-point...
t = 8;                                                                   % Prediction time-point
time_slope_px_pred = time_slope_px(strcmpi(time_slope_px.TP_diff, tp_diff(t-1)), :);	% PX change rate
px_pred = dose_time_px(strcmpi(dose_time_px.Time_point, time_pts(t)), :);           % PX data for pred tp
px_pred0 = dose_time_px(strcmpi(dose_time_px.Time_point, time_pts(t-1)), :);       % PX data for pred tp-1

resp_pred = dose_time_resp;                               % Response data
noise_types = {'None', 'Uniform', 'Gaussian'};
fprintf('What kind of noise do you want in response data?\n')
noise = noise_types{input('Enter 0, 1 or 2 for either "None", "Uniform" or "Gaussian".\nAns = ') + 1};
rng(0);         n_para = 0.05;                                  % Additive noise generator seed & parameter
switch lower(noise)
    case 'none'
        noise_resp = 0;
    case 'uniform'
        noise_resp = n_para * (- 1 + 2*rand(T * D * m, 1));
    case 'gaussian'
        noise_resp = n_para * (0 + 0.5*randn(T * D * m, 1));
end
resp_pred.Response = resp_pred.Response + noise_resp;   % Noisy resp

% m-fold crossvalidation...
train_idx = cell(m, 3);      test_idx = cell(m, 3);     Ytest = zeros(D, m);
Ypred = struct('HM', zeros(D, m), 'RF', zeros(D, m));
for i = 1 : m
    test_idx{i, 1} = (px_pred.Indiv == i);                         train_idx{i, 1} = ~test_idx{i, 1};
    test_idx{i, 2} = (time_slope_px_pred.Indiv == i);       train_idx{i, 2} = ~test_idx{i, 2};
    test_idx{i, 3} = (resp_pred.Indiv == i);                      train_idx{i, 3} = ~test_idx{i, 3};
    xx_train = px_pred{train_idx{i, 1}, 3:end};                xx0_train = px_pred0{train_idx{i, 1}, 3:end};
    xx_test = px_pred{test_idx{i, 1}, 3:end};                   xx0_test = px_pred0{test_idx{i, 1}, 3:end};
    % beta_train = time_slope_px_pred{train_idx{i, 2}, 3:end};
    % beta_test = time_slope_px_pred{test_idx{i, 2}, 3:end};
    % beta_train = [repmat(doses, [m-1, 1]), (xx_train(:, 2:end) - xx0_train(:, 2:end))];
    % beta_test = [doses, (xx_test(:, 2:end) - xx0_test(:, 2:end))];
    yy_train = reshape(resp_pred.Response(train_idx{i, 3}), [T, (m - 1)*D])';
    yy_test = reshape(resp_pred.Response(test_idx{i, 3}), [T, D])';
    
    % Form 3D ((m x D) x P x T) time series matrices...
    XX_train = cat(3, xx0_train(:, 2:end), xx_train(:, 2:end));
    XX_test = cat(3, xx0_test(:, 2:end), xx_test(:, 2:end));
    
    % Prediction...   
    tstep_ = 0;                                                       % One-step prediction => t_ = 0, t = 1
    n_tree = 100;     seed = 1;                                 % Use 'shuffle' for random outcomes
    ParameterModels = RecursiveHybridModel(doses, XX_train, yy_train, tstep_, n_tree, seed);        % HM
    yy_pred1 = ModelPredict(ParameterModels, XX_test, yy_test(:, t-1));
    rng(seed);          RFind = TreeBagger(n_tree, xx_train, yy_train(:, t), 'method', 'regression');      % RF
    yy_pred2 = predict(RFind, xx_test);
    
    Ytest(:, i) = yy_test(:, t);                                  % Validation matrix for pred tp
    Ypred.HM(:, i) = yy_pred1;                 Ypred.RF(:, i) = yy_pred2;
end

% % Error result table...
fprintf('Results summary: \n')
n_level = [min(noise_resp), max(noise_resp), mean(noise_resp), var(noise_resp)];
fprintf('Noise level: %s, n_min = %0.4f, n_max = %0.4f, n_mean = %0.4f, n_var = %0.4f\n', noise, n_level)
mspe = [mean((Ytest - Ypred.HM).^2, 1); mean((Ytest - Ypred.RF).^2, 1)]';     mspe = [mspe; mean(mspe, 1)];
mspe_table = array2table(round(mspe, 6), 'variablenames', {'MSPE_Hybrid', 'MSPE_Individual'},...
                                                    'rownames', [cellstr([repmat('S', [m, 1]), num2str((1:m)')]); {'Mean'}]);
disp(mspe_table)
fprintf('MSPE improvement = %0.4f (Hybrid / Individual)\n',...
                    mspe_table.MSPE_Hybrid(end) / mspe_table.MSPE_Individual(end))


%%
% Figures...
h = [figure('position', [50, 200, 900, 600]); figure('position', [950, 200, 900, 600])]';
mbox = setdiff(1:8, 4);      tym = str2num(cell2mat(split(time_pts, 't')));        y3D = zeros(T, D, m);    %#ok
for i = 1 : m
    y3D(:, :, i) = reshape(resp_pred{(resp_pred.Indiv == i), 4}, [T, D]);
    figure(h(1).Number),      subplot(2, 4, mbox(i))
    surf(doses, tym, y3D(:, :, i)),     view(-45,30),       colormap summer
    xlabel('Dose', 'FontName', 'Book Antiqua', 'FontWeight', 'bold', 'FontSize', 10, 'Rotation', 45)
    ylabel('Time', 'FontName', 'Book Antiqua', 'FontWeight', 'bold', 'FontSize', 10, 'Rotation', -30)
    zlabel('Response', 'FontName', 'Book Antiqua', 'FontWeight', 'bold', 'FontSize', 10)
    title(['Subject ', num2str(i)],...
        'FontName', 'Book Antiqua', 'FontWeight', 'bold', 'FontSize', 11, 'Color', [0, 0, 0])
    xticks(doses(1:2:end)),      yticks(tym(1:2:end));        zticks(0:0.2:1)
    xticklabels(strtrim(cellstr(num2str(round(doses(1:2:end), 2)))))
    yticklabels(cellstr(num2str(tym(1:2:end))))
    
    figure(h(2).Number),      subplot(2, 4, mbox(i))
    plot(doses, Ytest(:, i), 'gs--', doses, Ypred.HM(:, i), 'rs--', doses, Ypred.RF(:, i), 'bs--')
    xlabel('Dose ({\mu}M)', 'FontName', 'Book Antiqua', 'FontWeight', 'bold', 'FontSize', 11)
    ylabel('Response', 'FontName', 'Book Antiqua', 'FontWeight', 'bold', 'FontSize', 11)
    title(['Subject ', num2str(i)],...
        'FontName', 'Book Antiqua', 'FontWeight', 'bold', 'FontSize', 12, 'Color', [0, 0, 0])
    xticks(doses),      yticks(0:0.2:1)
    xticklabels(strtrim(cellstr(num2str(round(doses, 2)))))
end
ht = [ ];
titlelist = {'\bf3D Dose-Time Drug Response Surface'; '\bfDose-Response Curve Prediction'};
switch lower(noise)
    case 'uniform'
        addtitlepart = repmat({'\color[rgb]{0.5, 0, 0}Noise \sim {\itU}(\delta, \delta)'}, [2, 1]);
        titlelist = join([titlelist, addtitlepart], ', ');
    case 'gaussian'
        addtitlepart = repmat({'\color[rgb]{0.5, 0, 0}Noise \sim {\itN}(0, \delta^2/4)'}, [2, 1]);
        titlelist = join([titlelist, addtitlepart], ', ');
end
figure(h(1).Number),      ht{1} = suptitle(titlelist{1});
                                     set(ht{1}, 'FontName', 'Book Antiqua', 'FontWeight', 'bold', 'FontSize', 14)
figure(h(2).Number),     subplot(2,4,3),       legend({'Actual Data', 'Hybrid Model', 'Individual Model'},... 
                                                                                'Location', 'northeastoutside', 'FontName', 'Book Antiqua',...
                                                                                'FontWeight', 'bold', 'FontSize', 12)
                                     ht{2} = suptitle(titlelist{2});
                                     set(ht{2}, 'FontName', 'Book Antiqua', 'FontWeight', 'bold', 'FontSize', 14)

