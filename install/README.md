
# Install Support
This directory contains installation support files for the repository.

The `install/scripts/` directory contains installation-specific entrypoints and helper installers.

The `install/envs/` directory contains Conda/Mamba environment definitions used by `install/scripts/install_envs.sh`.

Creating `env_protocol` also installs this repository as an editable package with the selected manager. Reusing an existing `env_protocol` refreshes that editable installation. The explicit `--if_exists update` mode reconciles a YAML-backed environment to its YAML, installing declared packages and changing an installed version wherever the YAML declares a different one, which may mean a downgrade, and then performs the same refresh; packages the YAML no longer declares are left in place rather than pruned. Repeatable `--update_package` values bound that transaction to exact specifications already present in the YAML, and imply `--if_exists update`. Other environments are not package-install targets. Dry-run mode prints the environment and applicable package commands without executing them. Before/after transaction capture is a validation responsibility; the installer does not automatically preserve an environment-delta report.
