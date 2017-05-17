rm bardata
touch bardata

tot=0
out=""
seq=""
para=""

for tests in {1..4}
do
	case $tests in
		1)	printf '###\tseq\t1td\t2td\t3td\t4td\nSinglegrid100\t' >> bardata
		seq='./sequential 100 1500'
		para='./parallel 100 150000'
		;;
		2)	printf '\nSinglegrid200\t' >> bardata
		seq='./sequential 200 37500'
		para='./parallel 200 37500'
		;;
		3)	printf '\nMultigrid12\t' >> bardata
		seq='./mSeq 12 10400000'
		para='./mPara 12 10400000'
		;;
		4)	printf '\nMultigrid24\t' >> bardata
		seq='./mSeq 24 2630000'
		para='./mPara 24 2630000'
		;;
	esac

	echo -e "**** TEST $tests ****"
	echo "sequential:"
	tot=0
	for i in {1..5}
	do
		out=$($seq|grep Jacobi|cut -c "28-36")
		tot=$(echo $out + $tot | bc)
		echo $out
	done
	tot=$(bc <<< "scale=6; $tot/5")
	printf "%s\t""$tot"  >> bardata
	echo ""
	for n in {1..4}
	do
		echo "parallel ($n threads):"
		tot=0
		for i in {1..5}
		do
			out=$($para $n|grep Jacobi|cut -c "28-36")
			tot=$(echo $out + $tot | bc)
			echo $out
		done
		tot=$(bc <<< "scale=6; $tot/5")
		printf "%s\t""$tot"  >> bardata
		echo ""
	done
done
