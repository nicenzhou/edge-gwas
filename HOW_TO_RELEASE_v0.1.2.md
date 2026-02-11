# How to push edge-gwas v0.1.2 to GitHub and publish a release

Your repo already has the v0.1.2 commit and tag locally. Follow these steps to push to `main` and create the GitHub release.

---

## 1. (Optional) Commit any extra changes

You have uncommitted changes (new tests, PDF helper, etc.). To include them in v0.1.2:

```bash
cd /Users/nicen/edge/edge-gwas

# Stage what you want
git add .gitignore
git add tests/test_all_functions.py tests/test_all_functions_quick.py tests/test_report_pdf.py
git add edge_gwas/io_handlers.py edge_gwas/utils.py
git add RELEASE_v0.1.2.md HOW_TO_RELEASE_v0.1.2.md

# If you removed the old PDFs on purpose:
git add -u tests/test_run_essential.pdf tests/test_run_extend_1.pdf tests/test_run_extend_2.pdf tests/test_run_extend_3.pdf

git commit -m "Add full/quick test scripts with PDF reports; minor fixes"
```

If you want the tag v0.1.2 to point to this new commit, move the tag and force-push it (see step 3).

If you prefer to **keep the current tag** on the existing “Release v0.1.2” commit, skip the commit above and go to step 2.

---

## 2. Push `main` to GitHub

```bash
cd /Users/nicen/edge/edge-gwas
git push origin main
```

If Git asks for credentials, use your GitHub username and a **Personal Access Token** (not your password).  
Using SSH instead:

```bash
git remote set-url origin git@github.com:nicenzhou/edge-gwas.git
git push origin main
```

---

## 3. Push the tag `v0.1.2`

```bash
git push origin v0.1.2
```

That makes the tag visible on GitHub and allows you to create a release from it.

---

## 4. Create the GitHub Release (v0.1.2)

### Option A – GitHub website

1. Open: **https://github.com/nicenzhou/edge-gwas/releases**
2. Click **“Draft a new release”**.
3. **Choose a tag:** select **v0.1.2** (create from existing tag).
4. **Release title:** `v0.1.2`
5. **Description:** use the changelog for v0.1.2, for example:

```markdown
## v0.1.2 (2026-02-10)

### Bug fixes and improvements

- **io_handlers:** Added public API alias `format_gwas_output` so `save_for_locuszoom()` and package exports work.
- **utils:** Added `calculate_genomic_inflation(pvals)`; `create_summary_report()` now works.
- **core / io_handlers:** Replaced bare `except:` with `Exception` or `LinAlgError` for clearer errors.
- **setup.py:** Paths for requirements/README fixed; author/url set to Jiayan Zhou and nicenzhou/edge-gwas.
- Changelog and .gitignore added.

Install: `pip install git+https://github.com/nicenzhou/edge-gwas@v0.1.2`
```

6. Optionally attach **Source code (zip)** or a built wheel/sdist.
7. Click **“Publish release”**.

### Option B – GitHub CLI (`gh`)

If you have [GitHub CLI](https://cli.github.com/) installed and logged in:

```bash
cd /Users/nicen/edge/edge-gwas
gh release create v0.1.2 \
  --title "v0.1.2" \
  --notes "## v0.1.2 (2026-02-10)

### Bug fixes and improvements

- io_handlers: Added format_gwas_output alias; create_summary_report uses calculate_genomic_inflation.
- utils: Added calculate_genomic_inflation(pvals).
- core/io_handlers: Replaced bare except with Exception/LinAlgError.
- setup.py: Paths for requirements/README; author/url to Jiayan Zhou and nicenzhou/edge-gwas.

Install: pip install git+https://github.com/nicenzhou/edge-gwas@v0.1.2"
```

---

## 5. Verify

- **Main:** https://github.com/nicenzhou/edge-gwas (branch `main` shows latest commit).
- **Release:** https://github.com/nicenzhou/edge-gwas/releases (v0.1.2 listed).

Users can install this version with:

```bash
pip install git+https://github.com/nicenzhou/edge-gwas@v0.1.2
```

---

## Quick copy-paste (minimal: push existing v0.1.2)

If you only want to push the current commit and tag without new commits:

```bash
cd /Users/nicen/edge/edge-gwas
git push origin main
git push origin v0.1.2
```

Then create the release on the GitHub website (step 4, Option A).
