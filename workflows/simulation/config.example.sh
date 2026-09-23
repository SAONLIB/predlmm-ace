#!/usr/bin/env bash
# Copy to config.local.sh and edit for your cluster.

export OUT_ROOT="/path/to/predlmm/output"
export PYTHON_EXE="/path/to/python"
export GCTA_EXE="/path/to/gcta64"

export M=16000
export N_PCS=5
export BASE_SEED=20260921
export THREADS=8

export SLURM_ACCOUNT=""
export SLURM_PARTITION=""
export SLURM_CPUS=8
export SLURM_MEM="64G"
export SLURM_TIME="24:00:00"
export MAX_SIMULTANEOUS=30
