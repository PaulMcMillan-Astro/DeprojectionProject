samples="3 4 5"
bins=$(echo `python -c 'import numpy; print(str(numpy.logspace(-11, -9, 5))[1:-1])'`)


for sample in $samples; do
for bin in $bins; do
	printf '\n\n%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
	echo -e "Starting work on\nSample:                        Alpha:\nsynthetic_fivegauss_10e"$sample".vot   "$bin
	sed -i "26 c"$bin vars.ini
	sed -i "29 c synthetic/synthetic_fivegauss_10e"$sample".vot" vars.ini
	python synthetic_exe.py -i vars.ini -a 1
done
done

