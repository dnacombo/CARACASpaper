function [rej, failed_criteria] = CARACAS_rethresh(rej, CARACAS_struct, cfg_CARACAS)

% Initialize failed_criteria structure
failed_criteria = struct();

for i_c = 1:size(rej,2)
    NotCardiac = 0;
    
    % Initialize criteria failure flags for this component
    failed_criteria(i_c).sk = false;
    failed_criteria(i_c).ku = false;
    failed_criteria(i_c).PQ = false;
    failed_criteria(i_c).RR = false;
    failed_criteria(i_c).Rampl = false;
    failed_criteria(i_c).bpm = false;
    
    % Skewness
    if CARACAS_struct.meas(i_c).sk < cfg_CARACAS.thresh_sk
        NotCardiac = 1;
        failed_criteria(i_c).sk = true;
    end

    % Kurtosis
    if CARACAS_struct.meas(i_c).ku < cfg_CARACAS.thresh_ku
        NotCardiac = 1;
        failed_criteria(i_c).ku = true;
    end

    % RR
    if CARACAS_struct.meas(i_c).RR > cfg_CARACAS.thresh_RR
        NotCardiac = 1;
        failed_criteria(i_c).RR = true;
    end

    % Rampl
    if CARACAS_struct.meas(i_c).Rampl > cfg_CARACAS.thresh_Rampl
        NotCardiac = 1;
        failed_criteria(i_c).Rampl = true;
    end

    % bpm
    if CARACAS_struct.meas(i_c).bpm < cfg_CARACAS.thresh_bpm(1) || CARACAS_struct.meas(i_c).bpm > cfg_CARACAS.thresh_bpm(2)
        NotCardiac = 1;
        failed_criteria(i_c).bpm = true;
    end
    
    rej(1,i_c) = ~ NotCardiac;
end
end