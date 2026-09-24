# Unused env specs (audit D10, 2026-09-23)

No module, config or script references these. They are kept for reference only and are not
part of the pipeline or the paper's methods. If MaAsLin / ANCOM-BC become part of the analysis,
add their scripts and pin + lock the env like the others (conda_envs/lock/).

- ancombc_env.yml, maaslin_env.yml: unpinned; maaslin's header prescribes runtime installs
- krakentools_env.yml: empty (`dependencies: []`)
