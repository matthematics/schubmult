#!/bin/bash
# Run the C harness in parallel shards. Usage: run_c.sh harness_binary n outfile [workers]
set -e
BIN=$1; N=$2; OUT=$3; W=${4:-$(( $(nproc) - 2 ))}; MODE=${5:-}
TOTAL=1; for ((i=2;i<=N;i++)); do TOTAL=$((TOTAL*i)); done
CHUNK=$(( TOTAL / (W * 16) )); [ "$CHUNK" -lt 1 ] && CHUNK=1
rm -f "$OUT".part*
start=$(date +%s.%N)
for ((s=0; s<TOTAL; s+=CHUNK)); do
  e=$((s+CHUNK)); [ $e -gt $TOTAL ] && e=$TOTAL
  printf "%s %d %d %d %s.part%08d %s\n" "$BIN" "$N" "$s" "$e" "$OUT" "$s" "$MODE"
done | xargs -P "$W" -L 1 bash -c '"$0" "$1" "$2" "$3" "$4" $5'
cat "$OUT".part* > "$OUT"; rm -f "$OUT".part*
end=$(date +%s.%N)
echo "$BIN S_$N: $(( TOTAL * (TOTAL + 1) / 2 )) products, wall $(echo "$end - $start" | bc)s with $W workers -> $OUT"
