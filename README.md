##Finite difference method (explicit) by using CUDA - heat equation##

###Introduction###
This repo contains an implementation of finite difference method for heat equation with CUDA(Compute Unified Device Architecture).

###Environment###
- GPU : NVIDIA GeForce GTX 650 @ 1.072GHZ GDDR5 1GB
- [CUDA toolkit 7.5](https://developer.nvidia.com/cuda-toolkit)
- CPU : Intel(R) Core(TM) i5-6400 @ 2.7GHZ 
- RAM : DDR3L 16GB PC3-12800
- Microsoft Visual Studio Community 2013

###Result###
- `TESLA C2070`: faster than a CPU approximately 60 ~ 90 times 
- `NVIDIA GeForce GTX 650`: faster than a CPU approximately 10~20 times 

###Note###
- Current version is not optimized such as `shared memory`.
- If you're interested in my works, please visit my [homepage](https://sites.google.com/site/yoomh1989/).