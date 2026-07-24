
# Install Support
This directory contains installation support files for the repository.

The `install/scripts/` directory contains installation-specific entrypoints and helper installers.

The `install/envs/` directory contains Conda/Mamba environment definitions used by `install/scripts/install_envs.sh`.

Creating `env_protocol` also installs this repository as an editable package with the selected manager. Reusing an existing `env_protocol` refreshes that editable installation. The explicit `--if_exists update` mode reconciles declared YAML dependencies with `--freeze-installed`, without pruning or opportunistically updating installed packages, and then performs the same refresh. Repeatable `--update_package` values are the safer deliberate implementation for limiting a reviewed incremental transaction to exact specifications already present in the YAML. Other environments are not package-install targets. Dry-run mode prints the environment and applicable package commands without executing them. Before/after transaction capture is a validation responsibility; the installer does not automatically preserve an environment-delta report.
