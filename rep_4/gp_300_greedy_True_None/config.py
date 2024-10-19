cycle_dir: ''
diverse: true
metadata: Small test for active learning.
model_config:
  features:
    feature_type: fingerprint
    params:
      feature_column: Smiles
      fingerprint_radius: 3
      fingerprint_size: 2048
  hyperparameters: {}
  model_type: gp
  targets:
    feature_type: number
    params:
      feature_column: cnnaffinity
  tuning_hyperparameters:
  - {}
selection_config:
  hyperparameters: {}
  num_elements: 300
  selection_columns:
  - cnnaffinity
  - Smiles
  - h
  - fid
  selection_type: greedy
  tuning_hyperparameters: null
training_pool: ''
virtual_library: manual_init_h6_rgroups_linkers100_scorable.csv
time,1802.9842507839203,s