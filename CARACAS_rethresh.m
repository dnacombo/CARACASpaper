function withrej = CARACAS_rethresh(withrej, CARACAS_struct, cfg_CARACAS)

for i_c = 1:size(withrej,2)
    NotCardiac = 0;
    % Skewness
    if CARACAS_struct.meas(i_c).sk < cfg_CARACAS.thresh_sk
        NotCardiac = 1;
    end


    % RR
    if CARACAS_struct.meas(i_c).RR > cfg_CARACAS.thresh_RR
        NotCardiac = 1;
    end

    % bpm
    if CARACAS_struct.meas(i_c).bpm < cfg_CARACAS.thresh_bpm(1) || CARACAS_struct.meas(i_c).bpm > cfg_CARACAS.thresh_bpm(2)
        NotCardiac = 1;
    end


    % Kurtosis
    if CARACAS_struct.meas(i_c).ku < cfg_CARACAS.thresh_ku(1) || CARACAS_struct.meas(i_c).ku > cfg_CARACAS.thresh_ku(2)
        NotCardiac = 1;
    end

    % Rampl
    if CARACAS_struct.meas(i_c).Rampl > cfg_CARACAS.thresh_Rampl
        NotCardiac = 1;
    end

    withrej(1,i_c) = ~ NotCardiac;
end
end