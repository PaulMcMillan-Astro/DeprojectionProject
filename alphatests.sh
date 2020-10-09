sample_list="red blue solar full"
logalpha_list=$(echo `python -c 'import numpy; print(str(numpy.logspace(-11,-8,10))[1:-1])'`)

for sample in $sample_list; do
for logalpha in $logalpha_list; do
	echo "Currently: $sample $logalpha"
	sed -i "26 c$logalpha" vars.ini
	sed -i "29 c$sample""_sample.table" vars.ini
	python Deproject_exe.py -i vars.ini -a 1
done
done

