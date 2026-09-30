# Working with this repository

Repository: https://github.com/m-a-karbasian/composite-heat-transfer

## Get a working copy

```powershell
git clone https://github.com/m-a-karbasian/composite-heat-transfer.git
Set-Location composite-heat-transfer
```

Open this folder in MATLAB and follow the quick start in README.md.

## Publish later changes

From a clone of this repository, review and commit your changes:

```powershell
git status
git diff
git add .
git diff --cached --stat
git commit -m "Describe your changes"
git push
```

Use your own configured Git author identity and authenticate with your GitHub account when prompted. Never put credentials into project files.

The initial project files were published through the GitHub connector. An earlier local folder may have an unfinished Git initialization; use a fresh clone in a different folder to obtain the published history without overwriting that folder.

Review VALIDATION.md before describing numerical results or claiming that execution checks have passed. No license has been selected.
