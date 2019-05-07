if [ $# -lt 4 ]
then
	echo "sub.sh <SCRIPT> <BINS> <CONFIG> <SPID> "
else
	script=$1
	bins=$2
	cfg=$3
	spid=$4
	cmd="python3 "$script" "$cfg" "$bins" "$spid
	bsub -M 2000 -R "select[mem>2000] rusage[mem=2000]" -J  pd_"$spid" -eo pd_"$spid".e -oo pd_"$spid".o -q basement  "${cmd}"
fi
