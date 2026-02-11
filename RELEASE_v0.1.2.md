# Release v0.1.2 – Steps to Publish

The code is tagged and committed locally. Do the following to publish on GitHub.

## 1. Push to GitHub

From the repo root:

```bash
cd /Users/nicen/edge/edge-gwas
git push origin main
git push origin v0.1.2
```

(Use SSH if you have it: `git remote set-url origin git@github.com:nicenzhou/edge-gwas.git` then push.)

## 2. Create the GitHub Release

**Option A – GitHub website**

1. Open https://github.com/nicenzhou/edge-gwas/releases
2. Click **“Draft a new release”**
3. **Choose a tag:** select `v0.1.2` (create from existing tag).
4. **Release title:** `v0.1.2`
5. **Description:** use the changelog below (or paste from `docs/changelog.rst`).
6. Attach the source distribution: drag and drop `dist/edge-gwas-0.1.2.tar.gz` (from this repo).
7. Click **“Publish release”**.

**Option B – GitHub CLI**

If you have `gh` installed:

```bash
cd /Users/nicen/edge/edge-gwas
gh release create v0.1.2 dist/edge-gwas-0.1.2.tar.gz \
  --title "v0.1.2" \
  --notes-file - << 'EOF'
## v0.1.2 (2026-02-10)

### Bug fixes and improvements

- **io_handlers:** Added public API alias `format_gwas_output` so `save_for_locuszoom()` and package exports work.
- **utils:** Added `calculate_genomic_inflation(pvals)`; `create_summary_report()` now works.
- **core / io_handlers:** Replaced bare `except:` with `Exception` or `LinAlgError` for clearer errors.
- **setup.py:** Paths for requirements/README fixed; author/url set to Jiayan Zhou and nicenzhou/edge-gwas.
- Changelog and .gitignore added.

Install: `pip install git+https://github.com/nicenzhou/edge-gwas@v0.1.2`
EOF
```

## 3. Install from the new release

After the release is published, users can install with:

```bash
pip install git+https://github.com/nicenzhou/edge-gwas@v0.1.2
```

Or from PyPI if you publish there later:

```bash
pip install edge-gwas==0.1.2
```

## Built artifact

- **Source distribution:** `dist/edge-gwas-0.1.2.tar.gz` (already built in this repo).
