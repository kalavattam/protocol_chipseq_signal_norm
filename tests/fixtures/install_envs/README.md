
# Environment-installer fixtures
This fixture root generates fake `conda` and `mamba` commands for deterministic `install_envs.sh` state-transition contracts. The fakes record package-manager calls, maintain a caller-owned environment list, and can fail creation, installation, or editable-package refresh on request.

<br />

## Files
- `make.sh`: Generates the fake package-manager toolchain.
- `tool/`: Ignored generated command stubs.

<br />

## Current and deferred test coverage
The repository contract proves `fail`, `reuse`, and `update` command ordering, including successful creation, targeted update, and package-manager failure handling. A separate `RUN_INSTALL_ENVS=1` local integration uses a disposable real Conda environment to confirm creation and update against `env_siqchip`.
