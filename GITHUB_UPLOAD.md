# Upload to GitHub

After reviewing the author and license fields:

```bash
unzip predlmm-ace-v0.1.0.zip
cd predlmm-ace
git init
git branch -M main
git add .
git commit -m "Initial PredLMM-ACE research release"
```

The canonical repository is:

<https://github.com/SAONLIB/predlmm-ace>

To push to it with GitHub CLI authentication:

```bash
git remote add origin https://github.com/SAONLIB/predlmm-ace.git
git push -u origin main
```

or add the remote shown by GitHub:

```bash
git remote add origin git@github.com:SAONLIB/predlmm-ace.git
git push -u origin main
```

Keep the repository private until the manuscript authors approve the license,
software authorship, repository description, and public code/data statement.
Before a public release, add the repository URL and manuscript citation to
`CITATION.cff`, then create a versioned GitHub release and archive it in Zenodo.
