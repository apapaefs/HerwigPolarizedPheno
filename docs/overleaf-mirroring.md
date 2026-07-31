# GitHub and Overleaf mirroring

The manuscript was originally imported from Overleaf at commit `e37d50d`.
The current history deliberately remains a linear descendant of that commit;
an earlier unrelated-history replacement caused Overleaf's GitHub link to be
lost and must not be repeated.

This repository keeps `main.tex`, the bibliography, style files, and figures
at its root so the same commit is valid both as the complete phenomenology
repository and as the Overleaf paper bundle.

The intended remotes are:

```text
origin             git@github.com:apapaefs/HerwigPolarizedPheno.git
overleaf-github    git@github.com:apapaefs/Phenomenological-Investigations-of-Polarized-Collisions-in-Herwig-7.git
```

Before publishing a paper update, verify ancestry and then push the identical
`main` commit to both GitHub repositories:

```bash
git merge-base --is-ancestor e37d50d main
git push origin main
git push overleaf-github main:main
git ls-remote origin refs/heads/main
git ls-remote overleaf-github refs/heads/main
```

The two reported hashes must be identical. Overleaf's GitHub synchronization
can then pull from its existing linked repository without a history rewrite.

