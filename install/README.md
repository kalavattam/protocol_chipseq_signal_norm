
# Install Support
This directory contains installation support files for the repository.

The `install/scripts/` directory contains installation-specific entrypoints and helper installers.

The `install/envs/` directory contains Conda/Mamba environment definitions used by `install/scripts/install_envs.sh`.

Mamba `>= 1.5` or Conda `>= 24.7` is required, along with Bash `>= 4.4`. `install_envs.sh` stops with an explanatory message on anything older.

Creating `env_protocol` also installs this repository as an editable package with the selected manager. Reusing an existing `env_protocol` refreshes that editable installation. The explicit `--if_exists update` mode reconciles a YAML-backed environment to its YAML, installing declared packages and changing an installed version wherever the YAML declares a different one, which may mean a downgrade, and then performs the same refresh; packages the YAML no longer declares are left in place rather than pruned. Repeatable `--update_package` values bound that transaction to exact specifications already present in the YAML, and imply `--if_exists update`. Other environments are not package-install targets. Dry-run mode prints the environment and applicable package commands without executing them. Before/after transaction capture is a validation responsibility; the installer does not automatically preserve an environment-delta report.

<br />

## Installing where Conda channels are mirrored
For some institutions, public package channels are proxied and direct access to `anaconda.org` is blocked. In that case, standard channel names declared in `install/envs/*.yml` (e.g., `conda-forge` and `bioconda`) do not resolve, and an installation that does not name the local mirror cannot succeed at all. `install_envs.sh` (and `install_envs_entrypoint.sh`) allow users to pass the mirror to the installer:
```bash
bash install/scripts/install_envs.sh \
    --env_nam env_protocol \
    --channels "https://<mirror-host>/conda-forge/,https://<mirror-host>/bioconda/" \
    --override_channels \
    --yes
```

`--channels` searches the supplied channels ahead of those the environment YAML and `.condarc` declares, keeping the declared ones as a fallback. `--override_channels` removes that fallback: the YAML’s channels and the package manager’s defaults are both excluded, so only the supplied channels are searched. It takes no value of its own and requires `--channels`. Both apply to creation and to `--if_exists update`.

It’s recommended to prefer these flags over configuring channels globally, since they are explicit, recorded in the command that produced the environment, and leave the rest of a user’s package-manager configuration untouched. The latter may matter when installing across systems with and without blocked access to `anaconda.org`.

There are two things the installer does on behalf of the caller (both are reported by `--dry_run`, and both are retained upon failure, making it possible to inspect runs afterwards):
1. The installer renders the resulting channel list into a temporary copy of the environment YAML and installs from that copy; the tracked YAML is never modified.
2. It also renders a temporary package-manager configuration that disables channel redirection for that command alone.

The second thing matters for the following reason: a package manager may map a channel *name* onto a mirror list of its own, and because a mirrored name can match a supplied URL’s final path segment, the solved packages from the mirror can still be downloaded from a non-specified host.

Global configuration remains a separate matter. At a proxied site, Conda use *outside* this repository still reaches `anaconda.org`, so a `channel_alias` or a `.condarc` may still be wanted for that. It is neither needed by nor consulted for this installer.
