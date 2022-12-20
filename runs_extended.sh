#samples="blue red"
bins=$(echo `python -c 'import numpy; print(str(numpy.arange(0,10))[1:-1])'`)

#for sample in $samples; do
echo "##################################################"
echo "Running EDR3 with Alpha = 4e-10"
echo "##################################################"
sed -i "5 c""100, 100, 70" vars.ini
sed -i "8 c""-150,-150,-80" vars.ini

for bin in $bins; do
	echo "Currently: red_cut_""$bin"
	sed -i "29 c"500pc/3sigma/red_cut_"$bin"".vot" vars.ini
	python Deproject_exe.py -i vars.ini -a 1
done

for bin in $bins; do
	echo "Currently: blue_cut_""$bin"
	sed -i "29 c"500pc/3sigma/blue_cut_"$bin"".vot" vars.ini
	python Deproject_exe.py -i vars.ini -a 1
done


echo "##################################################"
echo "Running EDR3 with Alpha = 4e-10 and extended box"
echo "##################################################"
sed -i "5 c""150, 150, 70" vars.ini
sed -i "8 c""-225,-200,-80" vars.ini

for bin in $bins; do
	echo "Currently: red_cut_""$bin"
	sed -i "29 c"500pc/3sigma/red_cut_"$bin"".vot" vars.ini
	python Deproject_exe.py -i vars.ini -a 1
done

for bin in $bins; do
	echo "Currently: blue_cut_""$bin"
	sed -i "29 c"500pc/3sigma/blue_cut_"$bin"".vot" vars.ini
	python Deproject_exe.py -i vars.ini -a 1
done
