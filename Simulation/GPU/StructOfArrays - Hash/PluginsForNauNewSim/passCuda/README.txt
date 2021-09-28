1 - Gerar o projeto normalmente com o cmake
2 - nauPassCuda -> Properties -> Linker -> Input -> Additional Dependencies Pegar no caminho do tipo C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.2\lib\x64\cudart_static.lib Copiar e colar como C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.2\lib\x64\cudart.lib
3 - nauPassCuda -> Build Dependencies -> Build Customizations -> check CUDA 11.2
4 - mySort.cu -> Properties -> item Type cuda/c++
5 - Build