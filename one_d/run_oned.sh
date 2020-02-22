cd ../../build
make 
cd codes/one_d/

if [ "$1" = "solution" ]; then
	./main 
	python plot_sol.py
elif [ "$1" = "error_conv" ]; then
	./conv
else
	echo "Missing command line argument"
fi  
