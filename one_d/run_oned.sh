cd build
make 

if [ "$1" = "solution" ]; then
	./main 20 
	python plot_sol.py
elif [ "$1" = "error_conv" ]; then
	./conv
else
	echo "Missing command line argument"
fi  
