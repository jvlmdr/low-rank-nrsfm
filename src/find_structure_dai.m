function structure = find_structure_dai(projections, cameras, K)
  F = size(projections, 3);
  P = size(projections, 2);

  W = projections_to_matrix(projections);
  Rsh = full(block_diagonal_cameras(cameras));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Start of excerpt from NRSFM_BMM.m
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  S_PI = pinv(Rsh)*W;                       % Rank-3K S
  Shat_PI = S_to_Shat(S_PI,K);              % Transform S_PI to satisfy the K basis constraint

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ...
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  S_BMM = shape_recovery_fpca_s_sharp(W,Rsh,Shat_PI,K);

  Shat_BMM = S_to_Shat(S_BMM,K);                   % Transform to K basis form

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % End of excerpt from NRSFM_BMM.m
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  structure = structure_from_matrix(Shat_BMM);
end
