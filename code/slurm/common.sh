# 被 submit_*.sh / resubmit_failed.sh source：节点选择的公共设置
#   EXCLUDE   不使用的节点，逗号分隔（默认 n013；设 EXCLUDE= 为空则不排除）
#   N         数组作业数 = 同时占用的节点数。默认 = 分区里除 EXCLUDE 外所有可用节点（idle/mixed/allocated，排除 down/drain）
#   MAXN      N 的上限（默认 32）
#   EXTRA     其它 sbatch 参数，例如 EXTRA="--partition=xxx --account=yyy"
EXCLUDE=${EXCLUDE-n013}
EXTRA=${EXTRA:-}
[ -n "$EXCLUDE" ] && EXTRA="$EXTRA --exclude=$EXCLUDE"
CPUS=${CPUS:-512}
if [ -z "${N:-}" ]; then
    _part=$(echo "$EXTRA" | { grep -o -- '--partition=[^ ]*' || true; } | head -1 | cut -d= -f2)
    _nodes=""
    if command -v sinfo >/dev/null 2>&1; then
        _nodes=$(sinfo -h -N -t idle,mixed,allocated,completing ${_part:+-p "$_part"} -o "%N" 2>/dev/null | sort -u || true)
    fi
    if [ -n "$EXCLUDE" ] && [ -n "$_nodes" ]; then
        _nodes=$(echo "$_nodes" | { grep -v -x -F -f <(echo "$EXCLUDE" | tr ',' '\n') || true; })
    fi
    N=$(echo "$_nodes" | sed '/^$/d' | wc -l)
    [ "$N" -lt 1 ] && N=4
    [ "$N" -gt "${MAXN:-32}" ] && N=${MAXN:-32}
fi
echo "nodes: N=$N  cpus/task=$CPUS  exclude=${EXCLUDE:-none}  extra='$EXTRA'"
