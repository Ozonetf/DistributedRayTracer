# DistributedRayTracer
A RayTracer Engin project uisng hybrid architecture, MPI + PThread 

build:  mpicc *.cpp -lstdc++ -lpthread -lm \n
run:    mpirun -np "num processes" -hostfile ./cluster ./a.out

credits:
http://cosinekitty.com/raytrace/ \n
https://lodev.org/lodepng/
