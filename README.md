Eel
=====

"Fork Eel" to create a new MOOSE-based application.

For more information see: [http://mooseframework.org/create-an-app/](http://mooseframework.org/create-an-app/)


How to run: 

conda activate moose

make -j4

cd problems

../eel-opt -i input.i #replace input.i by whatever input file


cd ~/projects/babbler

mpiexec -n 4 ./babbler-opt -i test/tests/kernels/simple_diffusion/simple_diffusion.i

https://mooseframework.inl.gov/getting_started/examples_and_tutorials/tutorial01_app_development/step07_parallel.html
