% clc;    clear;      close all
clc;    clear

% Path & filenames...
PATH = 'Data/';
FILENAME1 = 'HMS_LINCS_RPPA_Data_Normalized_(Innoscan_Mapix)_SRD_Sep_21_azd.xlsx';
FILENAME2 = 'Copy of HMS_LINCS_Viability_Data_Normalized_SRD_Sep_21-2.xlsx';
FILENAME3 = 'HMS_LINCS_RPPA_Data_Extrapolated_All_72_hr_SRD_July_25.xlsx';

% Read data as Tables...
rppa_data = readtable([PATH, FILENAME1]);      rppa_data_header = rppa_data.Properties.VariableNames';
apop_data = readtable([PATH, FILENAME2]);     apop_data_header = apop_data.Properties.VariableNames';

drug_list = intersect(unique(rppa_data.SmallMoleculeName), unique(apop_data.SmallMoleculeName));
cl_list = intersect(unique(rppa_data.CellLineName), unique(apop_data.CellLineName));
dose_list = intersect(unique(rppa_data.DrugConcentration_uM_), unique(apop_data.DrugConcentration_uM_));
doses = dose_list(2:end);                                     % Dose #1 == 0 for DMSO control
tp_list_rppa = unique(rppa_data.TimePoint_hr_);     
tp_list_apop = unique(apop_data.TimePoint_hr_);
tp_list = intersect(tp_list_rppa, tp_list_apop);     % Common time-points betn rppa & apoptosis data
protein_list = rppa_data.Properties.VariableNames(5:end)';

% #parameters...
C = numel(drug_list);          P = numel(protein_list);
Tx = numel(tp_list_rppa);	Ty = numel(tp_list_apop);
T = numel(tp_list);              D = numel(doses);          m = numel(cl_list);

% Variable nomenclature...
delim = {' ', '-', '_'};
drug_names = [{'dmso'}; lower(drug_list(2:end))];          % 1st one is the DMSO control
for k = 2 : C
    drug_names{k} = cell2mat(join(split(drug_names{k}, delim), ''));
end

cl_names = lower(cl_list);
for i = 1 : m
    cl_names{i} = cell2mat(join(split(cl_names{i}, delim), ''));
end

% Drug-wise analysis...
% drug_idx = 2;                                                   % , AZ-628, PLX-4720, SB590885, Selumetinib, Vemurafenib
fprintf('Enter the desired drug index... \n')
fprintf('\t (AZ-628 = 2, PLX-4720 = 3, SB590885 = 4, Selumetinib = 5, Vemurafenib = 6)\n')
drug_idx = input('Index = ');
chosen_drug = drug_list{drug_idx};
chosen_drug_name = drug_names{drug_idx};
fprintf(1, 'Chosen drug for analysis = '),      fprintf(2, '%s\n', chosen_drug)

% Generate data variables for 3 time points...
var_names = cell(m, 2, Ty);
for t = 1 : Ty
    chosen_tp = tp_list_apop(t);                             % [24, 48, 72]
    
    % Drug-wise Tables (Creates dynamic variables... NOT advised in general!!!)...
    drug_table_names = {['rppa_data_', chosen_drug_name]; ['apop_data_', chosen_drug_name]};
    eval([drug_table_names{1}, ' = rppa_data(strcmpi(rppa_data.SmallMoleculeName, chosen_drug), :);']);
    eval([drug_table_names{2}, ' = apop_data(strcmpi(apop_data.SmallMoleculeName, chosen_drug), :);']);
    
    for i = 1 : m
        rppa_idx_cl_tp = eval(['intersect(find(strcmpi(', drug_table_names{1}, '.CellLineName, cl_list{i})), find(',...
                                                drug_table_names{1}, '.TimePoint_hr_ == ', num2str(chosen_tp), '));']);
        apop_idx_cl_tp = eval(['intersect(find(strcmpi(', drug_table_names{2}, '.CellLineName, cl_list{i})), find(',...
                                                drug_table_names{2}, '.TimePoint_hr_ == ', num2str(chosen_tp), '));']);
        
        % Create dynamic variables (NOT advised in general!!!)
        var_names(i, :, t) = {[cl_names{i}, '_t', num2str(chosen_tp), '_data'];
                                        [cl_names{i}, '_yy_t', num2str(chosen_tp), '_data']};
        eval([var_names{i, 1, t}, ' = ', drug_table_names{1}, '{rppa_idx_cl_tp, 5:end};']);
        eval([var_names{i, 2, t}, ' = ', drug_table_names{2}, '{apop_idx_cl_tp, 5:end};']);
    end
end

% Read extrapolated rppa data...
rppa_exp_data = readtable([PATH, FILENAME3]);
rppa_exp_data_header = rppa_exp_data.Properties.VariableNames';
drug_exp_table_name = ['rppa_exp_data_', chosen_drug_name];
eval([drug_exp_table_name, ' = rppa_exp_data(strcmpi(rppa_exp_data.SmallMoleculeName, chosen_drug), :);']);
for i = 1 : m
	rppa_exp_idx_cl_tp = eval(['intersect(find(strcmpi(', drug_exp_table_name,...
                                       '.CellLineName, cl_list{i})), find(', drug_exp_table_name, '.TimePoint_hr_ == 72));']);	
    eval([var_names{i, 1, 3}, ' = ', drug_exp_table_name, '{rppa_exp_idx_cl_tp, 5:end};']);
end

%%
% var_names(cl_name, rppa/apop, time_point)
 
% Model evaluation...
cl_test = {'K2'; 'MMAC-SF'; 'SKMEL28'};            [~, test_idx] = intersect(lower(cl_list), lower(cl_test));
cl_train = setdiff(cl_list, cl_test);                       [~, train_idx] = setdiff(1:m, test_idx);
m_test = numel(cl_test);                                       m_train = numel(cl_train);

xxtrain_t24 = eval([ '[', cell2mat(join(var_names(train_idx, 1, 1), '; ')), ']' ]);
xxtrain_t48 = eval([ '[', cell2mat(join(var_names(train_idx, 1, 2), '; ')), ']' ]);
xxtrain_t72 = eval([ '[', cell2mat(join(var_names(train_idx, 1, 3), '; ')), ']' ]);     % Extrapolated px, used for RF
xxtrain = cat(3, xxtrain_t24, xxtrain_t48);                                         % Concatenate in a 3D array

% disp([ '[', cell2mat(join(var_names(train_idx, 1, 1), '; ')), ']' ])          % Shows the evaluated command line

% Time-series data (Ensemble of response data for all 7 training CLs for a particular tp)...
yytrain_t24 = eval([ '[', cell2mat(join(var_names(train_idx, 2, 1), '; ')), ']' ]);
yytrain_t48 = eval([ '[', cell2mat(join(var_names(train_idx, 2, 2), '; ')), ']' ]);
yytrain_t72 = eval([ '[', cell2mat(join(var_names(train_idx, 2, 3), '; ')), ']' ]);
yytrain = [yytrain_t24, yytrain_t48, yytrain_t72];          % (m_train x D) x T

%% 
% Prediction...
mD_test = m_test * D;
xxtest_t24 = eval([ '[', cell2mat(join(var_names(test_idx, 1, 1), '; ')), ']' ]);
xxtest_t48 = eval([ '[', cell2mat(join(var_names(test_idx, 1, 2), '; ')), ']' ]);
xxtest_t72 = eval([ '[', cell2mat(join(var_names(test_idx, 1, 3), '; ')), ']' ]);     % Extrapolated px, used for RF
xxtest = cat(3, xxtest_t24, xxtest_t48);                                         % Concatenate in a 3D array

yytest_t24 = eval([ '[', cell2mat(join(var_names(test_idx, 2, 1), '; ')), ']' ]);
yytest_t48 = eval([ '[', cell2mat(join(var_names(test_idx, 2, 2), '; ')), ']' ]);
yytest_t72 = eval([ '[', cell2mat(join(var_names(test_idx, 2, 3), '; ')), ']' ]);
% yytest = [yytest_t24, yytest_t48, yytest_t72];          % (m_test x D) x T

tp_pred = [Ty-1, Ty, Ty];                                     % Prediction tp
yytest_names = reshape(var_names(test_idx, 2, :), [m_test, Ty]);
kklin = {(tp_pred - 2) * m_test + (1:m_test); (tp_pred - 1) * m_test + (1:m_test)};            % [t_; t] linear indices
yytest = eval(['[', cell2mat(join(yytest_names(kklin{2}), ', ')), ']']);            % D x m_test

tic
% Recursive Hybrid model...
tstep_ = 0;                                                           % One-step prediction => t_ = 0, t = 1
n_tree = 100;     seed = 1;                                     % Use 'shuffle' for random outcomes
ParameterModels = RecursiveHybridModel(doses, xxtrain, yytrain, tstep_, n_tree, seed);

yy0test = eval(['[', cell2mat(join(yytest_names(kklin{1}), '; ')), ']']);          % (m_test x D) x 1
yypred.HM = ModelPredict(ParameterModels, xxtest, yy0test);
toc

tic
% Individual RF models...
xxtrainRF_t72 = [repmat(doses, m_train, 1), xxtrain_t72];
xxtrainRF_t48 = [repmat(doses, m_train, 1), xxtrain_t48];
rng(seed);          RF72 = TreeBagger(n_tree, xxtrainRF_t72, yytrain_t72, 'method', 'regression');
rng(seed);          RF48 = TreeBagger(n_tree, xxtrainRF_t48, yytrain_t48, 'method', 'regression');      % For K2

xxtest_names = reshape(var_names(test_idx, 1, :), [m_test, Ty]);
RFind = cellfun(@eval, cellstr([repmat('RF', [m_test, 1]), num2str(tp_list_apop(tp_pred))]), 'uniformoutput', 0);
xxtestRF_tp = zeros(D, P+1, m_test);
for i = 1 : m_test
    xxtestRF_tp(:, :, i) = [doses, eval(cell2mat(xxtest_names(kklin{2}(i))))];
    yypred.RF(:, i) = predict(RFind{i}, xxtestRF_tp(:, :, i));
end
toc

% Error calculation...
mspe = [mean((yytest - yypred.HM).^2, 1)', mean((yytest - yypred.RF).^2, 1)'];      mspe = [mspe; mean(mspe, 1)];
mspe_table = array2table(round(mspe, 6), 'variablenames', {'MSPE_Hybrid', 'MSPE_Individual'},...
                                                                                            'rownames', [cl_test; {'Mean'}]);
disp(mspe_table)
fprintf('MSPE improvement = %0.4f (Hybrid / Individual)\n',...
                                                mspe_table.MSPE_Hybrid(end) / mspe_table.MSPE_Individual(end))


%% 
% Figures...
figure('position', [500, 200, 900, 700])
yrange = zeros(1, m_test);
for k = 1 : m_test
    subplot(m_test,1,k)
    plot(doses, yytest(:, k), 'gs--', doses, yypred.HM(:, k), 'rs--', doses, yypred.RF(:, k), 'bs--')
    xlabel('Dose ({\mu}M)', 'FontName', 'Book Antiqua', 'FontWeight', 'bold', 'FontSize', 11)
    ylabel({'Apoptosis'; 'Fraction'}, 'FontName', 'Book Antiqua', 'FontWeight', 'bold', 'FontSize', 11)
    title(['Cell line = ', cl_test{k}], 'FontName', 'Book Antiqua', 'FontWeight', 'bold', 'FontSize', 12)
	yrange(k) = round(max([yytest(:, k); yypred.HM(:, k); yypred.RF(:, k)]), 1) + ((k == 1)*0.1 + (k > 1)*0.05);
    axis([0, 3.2, 0, yrange(k)]),      grid on
    xticks(0:0.4:3.2),      yticks(linspace(0, yrange(k), 5)),      yticklabels(round(yticks, 2))
end
subplot(m_test,1,1),        legend({'Actual Data', 'Hybrid Model', 'Individual Model'},... 
                                                    'Location', 'north', 'Orientation', 'horizontal', 'FontName', 'Book Antiqua',...
                                                    'FontWeight', 'bold', 'FontSize', 12)
ht = suptitle(['Dose - Response Curve Prediction, Drug = \color[rgb]{0.5, 0, 0}', drug_list{drug_idx}]);
set(ht, 'FontName', 'Book Antiqua', 'FontWeight', 'bold', 'FontSize', 14)
