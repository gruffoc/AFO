for i in $(ls *.dat); do
	
	a=$(echo $i | awk -F "_" {'print $2'} | awk -F "." {'print $1'})
	c=$(cat $i)
	
	echo $a $c;

done

