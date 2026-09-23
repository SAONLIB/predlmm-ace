# Contributing

Please open an issue before changing the statistical model or output schema.
Changes should include a focused test and must preserve exact sample-order
validation. Run:

```bash
python -m pip install -e '.[test]'
pytest
bash workflows/simulation/smoke_test.sh
```

Do not commit controlled data, participant identifiers, local cluster paths,
large GRM/genotype files, outputs, or credentials.
