function [Cy, Cx11, Cx12, Cx21, Cx22, obj] = scICML_main(p11,p12,p21,p22,Cy, Cx11, Cx12, Cx21, Cx22, alpha,iter,nh)
obj = zeros(iter,1);
%tic
[tilde_p11, cluster_p11] = updateTildep_plus(p11, Cx11, Cy);
[tilde_p12, ~] = updateTildep_plus(p12, Cx12, Cy);
[tilde_p21, cluster_p21] = updateTildep_plus(p21, Cx21, Cy);
[tilde_p22, ~] = updateTildep_plus(p22, Cx22, Cy);
[ind_X11, ind_X21]=swap_label_plus(cluster_p11,cluster_p21,nh);
%toc

for t = 1:iter
  tic
  % upate the coefficient of the matching term D_{KL}(...||...)
  Coef = updateMatchCoef(cluster_p11, cluster_p21, ind_X11, ind_X21, alpha);
  
  % update Cxijbased on pij, tilde_pij, Coef;
  Cx11 = updateRowClustering_X(p11, tilde_p11, Cx11, Coef);
  Cx12 = updateRowClustering_pc(p12, tilde_p12, Cx12);
  Cx21 = updateRowClustering_X(p21, tilde_p21, Cx21, Coef);
  Cx22 = updateRowClustering_pc(p22, tilde_p22, Cx22);

  % update Cy
  [Cy,dist] = updateColClustering_Y(p11, p12, p21, p22, tilde_p11, tilde_p12, tilde_p21, tilde_p22, Cy, alpha, Coef);

  %show the ARI NMI RAGI;
  %printResult_lv(Cy,data_score, gnd);

  % update tilde_pijbased on pij, Cxij, Cy
  [tilde_p11, cluster_p11] = updateTildep_plus(p11, Cx11, Cy);
  [tilde_p12, ~] = updateTildep_plus(p12, Cx12, Cy);
  [tilde_p21, cluster_p21] = updateTildep_plus(p21, Cx21, Cy);
  [tilde_p22, ~] = updateTildep_plus(p22, Cx22, Cy);

  % match the cell types across the target data and the source data;
  [ind_X11, ind_X21]=swap_label_plus(cluster_p11,cluster_p21,nh);
  %matm(t,:) = [ind_X11,ind_X21]; %record the matching history;
  
  % record values of objective function in each update
  obj(t) = dist + Coef;
  toc
end