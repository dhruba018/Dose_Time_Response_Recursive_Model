clc;    clear;      close all

% Path & filenames...
PATH = 'Data/';
FILENAME1 = 'HMS_LINCS_RPPA_Data_Normalized_(Innoscan_Mapix)_SRD_Sep_21_azd.xlsx';

% Read data as table...
pData = readtable([PATH, FILENAME1]);       pDataHeader = pData.Properties.VariableNames';
timepts = unique(pData.TimePoint_hr_);          proteins = pDataHeader(5:end);
clines = unique(pData.CellLineName);
drug_list = unique(pData.SmallMoleculeName);
dose_list = unique(pData.DrugConcentration_uM_);

% Discarding DMSO control data...
drugs = drug_list(2:end);       
doses = dose_list(2:end);                % Discarding dose == 0

% 72-hr data matrix...
timept_exp = 72;
nCD = numel(clines) * numel(doses);
pDataDrug72 = cell(nCD, numel(pDataHeader), numel(drugs));
pData72 = [ ];
for k = 1:numel(drugs)
    pDataDrug = pData(strcmpi(pData.SmallMoleculeName, drugs{k}), :);
    pDataDrug72(:, 1, k) = reshape(repmat(clines', [numel(doses), 1]), [nCD, 1]);       % CL column
    pDataDrug72(:, 2, k) = repmat(drugs(k), [nCD, 1]);                                              % Drug column
    pDataDrug72(:, 3, k) = num2cell(repmat(doses, [numel(clines), 1]));                     % Dose column
    pDataDrug72(:, 4, k) = repmat(num2cell(timept_exp), [nCD, 1]);                          % TP column
    for i = 1:numel(clines)
        for j = 1:numel(doses)
            lineDoseData = pDataDrug(strcmpi(pDataDrug.CellLineName, clines(i)) &...
                                        pDataDrug.DrugConcentration_uM_ == doses(j), :);         % Data for dose - CL pairs
            xx = table2array(lineDoseData(:, 5:end));
            xx_exp = spline(timepts, xx', timept_exp)';                                                % Spline extrapolation
            pDataDrug72((i - 1)*numel(doses) + j, 5:end, k) = num2cell(xx_exp);
        end
    end
    
    % Construct 72-hr data table...
    pData72 = [pData72; cell2table(pDataDrug72(:, :, k), 'variablenames', pDataHeader)];        %#ok
    
end

% % Writing in file...
% FILENAME2 = 'HMS_LINCS_RPPA_Data_Extrapolated_All_72_hr_SRD2_July_25.xlsx';
% writetable(pData72, [PATH, FILENAME2], 'writevariablenames', 1, 'sheet', '72_hr_spline_extrapolation')
% RemoveExcelSheet([PATH, FILENAME2])            % Removes defalut excel sheets

