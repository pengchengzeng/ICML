function ac = NMIARI(indic, label)
%X=V_ini{i}; K=K0; kmeansFlag=options.kmeans;
%X = V_scAI.'; data_score= data_score; label= gnd;kmeansFlag=1;
[~, nmi_value, ~] = CalcMetrics(label, indic);
[ari_value,~,~,~] = RandIndex(label,indic);
fprintf('NMI:%0.4f\t\tARI:%0.4f\t\n', nmi_value, ari_value);