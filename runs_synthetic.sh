alphas=$(echo `python -c 'import numpy; print(str(numpy.logspace(-13, -11, 5)[:-1])[1:-1])'`)
IFS=', ' read -r -a alphas <<< "$alphas"
for i in {0..3}; do
	echo "Currently on: ${alphas[i]}" > logs/alpha_test$i.log
	sed -i "26 c"${alphas[i]} vars.ini
    tail vars.ini | head -n 2
	python Deproject_exe.py -i vars.ini -a 1 -b 0 -bs 0 >> logs/alpha_test$i.log &
    sleep 60
done

