#!/bin/bash
# check_colonization_progress.sh — snapshot progress of a colonization
# run_colonization.R job (single run or sweep) from its log file.
# Read-only — safe to run any time, no module load needed.
#
# Usage: sh scripts/diagnostics/check_colonization_progress.sh <site> <exp_tag>
#   e.g. sh scripts/diagnostics/check_colonization_progress.sh Maquipucuna pollination_success
# ─────────────────────────────────────────────────────────────────────────────
SITE=$1
TAG=$2
if [ -z "$SITE" ] || [ -z "$TAG" ]; then
  echo "Usage: $0 <site> <exp_tag>"
  exit 1
fi

LOG=$(ls -t logs/colonization_${SITE}_${TAG}_*.log 2>/dev/null | head -1)
if [ -z "$LOG" ]; then
  echo "No log found matching logs/colonization_${SITE}_${TAG}_*.log"
  exit 1
fi
echo "Log: $LOG"
echo

SWEEP_LINE=$(grep -E "Sweeping|Factorial sweep" "$LOG")
if [ -n "$SWEEP_LINE" ]; then
  echo "$SWEEP_LINE"
  # One-at-a-time: "N values ... with M reps"; factorial: "N combos x M reps"
  N_VALUES=$(echo "$SWEEP_LINE" | grep -oE '[0-9]+ values' | grep -oE '^[0-9]+')
  N_COMBOS=$(echo "$SWEEP_LINE" | grep -oE '[0-9]+ combos' | grep -oE '^[0-9]+')
  N_REPS=$(echo "$SWEEP_LINE" | grep -oE '[0-9]+ reps' | grep -oE '^[0-9]+')
  N_CORES=$(echo "$SWEEP_LINE" | grep -oE '[0-9]+ cores' | grep -oE '^[0-9]+')
  N_CORES=${N_CORES:-"?"}
  if [ -n "$N_COMBOS" ] && [ -n "$N_REPS" ]; then
    N_JOBS=$((N_COMBOS * N_REPS))
  elif [ -n "$N_VALUES" ] && [ -n "$N_REPS" ]; then
    N_JOBS=$((N_VALUES * N_REPS))
  else
    N_JOBS="?"
  fi
  DONE_STEPS=$(grep -c '\] t=' "$LOG")
  echo "Timestep log lines so far: $DONE_STEPS  (all $N_JOBS jobs interleaved in this file — only $N_CORES run concurrently)"
else
  DONE_STEPS=$(grep -c '\] t=' "$LOG")
  echo "Single (non-swept) run — timestep log lines so far: $DONE_STEPS"
fi

N_ERR=$(grep -c '^ERROR' "$LOG")
if [ "$N_ERR" -gt 0 ]; then
  echo
  echo "$N_ERR worker(s) reported errors:"
  grep '^ERROR' "$LOG"
fi

echo
echo "--- last 10 log lines ---"
tail -10 "$LOG"

echo
echo "--- your squeue ---"
squeue -u "$USER"
