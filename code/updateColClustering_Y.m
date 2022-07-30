function [Cz, dist] = updateColClustering_Y(p, q, p1, q1, tilde_p, tilde_q, tilde_p1, tilde_q1, Cz, lambda, Coef)
%p=p11;q=p12; p1=p21; q1=p22; tilde_p=tilde_p11; tilde_q=tilde_p12; tilde_p1=tilde_p21; tilde_q1=tilde_p22; Cz=Cy; lambda=1;
% compute p(X|z),p1(X|z) => p11, p12
pXz = p./repmat(sum(p,1), size(p,1), 1);
p1Xz = p1./repmat(sum(p1,1), size(p1,1), 1);

% compute q(Y|z),q1(Y|z) => p21,p22
qYz = q./repmat(sum(q,1), size(q,1), 1);
q1Yz = q1./repmat(sum(q1,1), size(q1,1), 1);

% compute tilde_p(X|tilde_z)
tilde_pXz = tilde_p./repmat(sum(tilde_p,1), size(tilde_p,1), 1);

for i = 1:size(tilde_p,1)
  tilde_pXtz(i,:) = accumarray(Cz, tilde_p(i,:)')';
end

% compute tilde_p1(X|tilde_z)
tilde_p1Xz = tilde_p1./repmat(sum(tilde_p1,1), size(tilde_p1,1), 1);

for i = 1:size(tilde_p1,1)
  tilde_p1Xtz(i,:) = accumarray(Cz, tilde_p1(i,:)')';
end

% compute tilde_q(Y|tilde_z)
tilde_qYz = tilde_q./repmat(sum(tilde_q,1), size(tilde_q,1), 1);

for i = 1:size(tilde_q,1)
  tilde_qYtz(i,:) = accumarray(Cz, tilde_q(i,:)')';
end

% compute tilde_q1(Y|tilde_z)
tilde_q1Yz = tilde_q1./repmat(sum(tilde_q1,1), size(tilde_q1,1), 1);

for i = 1:size(tilde_q1,1)
  tilde_q1Ytz(i,:) = accumarray(Cz, tilde_q1(i,:)')';
end

% find Cz minimizing objective function
pz = sum(p,1); qz = sum(q,1);
p1z = sum(p1,1); q1z = sum(q1,1);

for zc = 1:size(tilde_pXtz,2)
  for z = 1:size(pXz,2)
    temp_X11(zc, z) = pz(z) * KLDiv(pXz(:,z)', tilde_pXtz(:,zc)');
    temp_X12(zc, z) = qz(z) * KLDiv(qYz(:,z)', tilde_qYtz(:,zc)');
    temp_X21(zc, z) = p1z(z) * KLDiv(p1Xz(:,z)', tilde_p1Xtz(:,zc)');
    temp_X22(zc, z) = q1z(z) * KLDiv(q1Yz(:,z)', tilde_q1Ytz(:,zc)');    
  end
end

temp = temp_X11 + lambda * temp_X12 + temp_X21 + temp_X22 + Coef / size(p,2);

[mindist, Cz] = min(temp);
Cz = Cz';
dist = sum(mindist);

clearvars -except Cz dist
