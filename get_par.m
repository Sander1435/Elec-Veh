dialogPrms = get_param('EV_2spd_4AUB10_2023b/Multispeed gearbox', 'DialogParameters');
dialogPrmNames = fieldnames(dialogPrms);
for idx = 1:numel(dialogPrmNames)
    prmName = dialogPrmNames{idx};
    prmValue = get_param('EV_2spd_4AUB10_2023b/Multispeed gearbox', prmName);
    fprintf('Parameter Name: %s, Value: %s\n', prmName, prmValue);
end
