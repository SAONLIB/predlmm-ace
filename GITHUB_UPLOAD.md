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

Create an empty GitHub repository named `predlmm-ace`, then either use GitHub
CLI:

```bash
gh repo create predlmm-ace --private --source=. --remote=origin --push
```

or add the remote shown by GitHub:

```bash
git remote add origin git@github.com:YOUR_ACCOUNT/predlmm-ace.git
git push -u origin main
```

Keep the repository private until the manuscript authors approve the license,
software authorship, repository description, and public code/data statement.
Before a public release, add the repository URL and manuscript citation to
`CITATION.cff`, then create a versioned GitHub release and archive it in Zenodo.
