clc;    clear;      close all

pData = readtable('Data\HMS_LINCS_RPPA_Data_Normalized_(Innoscan_Mapix)_SRD_Sep_21_azd.xlsx');
drugs = unique(pData.SmallMoleculeName);
proteins = pData.Properties.VariableNames(5:end)';
pData_azd = pData(strcmpi(pData.SmallMoleculeName, 'AZ-628'), :);
doses = unique(pData_azd.DrugConcentration_uM_);
clines = unique(pData_azd.CellLineName);
timepts = unique(pData_azd.TimePoint_hr_);

pData_azd_72 = cell(numel(clines) * numel(doses), numel(proteins) + 4);
pData_azd_72(:, 1) = reshape(repmat(clines', [numel(doses), 1]), [numel(clines)*numel(doses), 1]);
pData_azd_72(:, 2) = repmat({'AZ-628'}, [numel(clines) * numel(doses), 1]);
pData_azd_72(:, 3) = num2cell(repmat(doses, [numel(clines), 1]));
pData_azd_72(:, 4) = repmat(num2cell(72), [numel(clines) * numel(doses), 1]);
for i = 1:numel(clines)
    for j = 1:numel(doses)
        lineDoseData = pData_azd(strcmpi(pData_azd.CellLineName, clines(i)) &... 
                                                    pData_azd.DrugConcentration_uM_ == doses(j), :);
        xx = table2array(lineDoseData(:, 5:end));
        xx_exp = spline(timepts, xx', 72)';
        pData_azd_72((i - 1)*numel(doses) + j, 5:end) = num2cell(xx_exp);
    end
end

pData_azd_72 = cell2table(pData_azd_72, 'VariableNames', pData_azd.Properties.VariableNames);
% writetable(pData_azd_72, 'Data\HMS_LINCS_RPPA_Data_Extrapolated_72_hr_SRD_June_17_azd.xlsx',...
%                     'WriteVariableNames', 1)




