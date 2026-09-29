#!/bin/bash
# batch_summary_mailer.sh — one-page batch summary, generic across experiments.
#
# Not itself an experiment script -- a reusable tail job. Submit with
# --dependency=afterany:<jobid1>:<jobid2>:... and --mail-type=END so SLURM's
# own completion notice for THIS job doubles as "the batch is done"; the
# actual summary content (job name/state/elapsed/exit code/output path per
# job in the batch) is composed here from sacct and mailed directly, since
# SLURM's own --mail-type mail has fixed content and can't carry it.
#
# Usage: sbatch --dependency=afterany:JID1:JID2:... --job-name=<experiment>_summary \
#          --mail-type=END --mail-user=s38leste@uni-bonn.de \
#          --output=logs/%x_%j.out --error=logs/%x_%j.err \
#          scripts/diagnostics/batch_summary_mailer.sh <BATCH_LABEL> JID1 JID2 ...
#
# Run from: /home/s38leste_hpc/seres/
#SBATCH --partition=intelsr_short
#SBATCH --account=ag_biob_scabral
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
set -uo pipefail
cd /home/s38leste_hpc/seres
BATCH_LABEL="$1"; shift
JIDS=("$@")

SUMMARY=$(mktemp)
{
  echo "Batch summary: ${BATCH_LABEL}"
  echo "Generated: $(date)"
  echo ""
  printf "%-12s %-30s %-12s %-10s %-10s %s\n" "JobID" "JobName" "State" "ExitCode" "Elapsed" "OutputPath"
  for jid in "${JIDS[@]}"; do
    line=$(sacct -j "$jid" -X --noheader --parseable2 --format=JobID,JobName,State,ExitCode,Elapsed)
    IFS='|' read -r id name state exitcode elapsed <<< "$line"
    outpath=$(scontrol show job "$jid" 2>/dev/null | grep -oP 'StdOut=\K\S+')
    [ -z "$outpath" ] && outpath="(job purged from scontrol -- see logs/*_${jid}.out)"
    printf "%-12s %-30s %-12s %-10s %-10s %s\n" "$id" "$name" "$state" "$exitcode" "$elapsed" "$outpath"
  done
} > "$SUMMARY"

cat "$SUMMARY"
mail -s "[seres] Batch summary: ${BATCH_LABEL}" s38leste@uni-bonn.de < "$SUMMARY"
rm -f "$SUMMARY"
