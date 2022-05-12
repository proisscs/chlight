This is the code used for the experiments of the publication "Light Contraction Hierarchies: Hierarchical Search Without Shortcuts"

Instructions for Ubuntu 20.04:

Step 1: clone the repository and navigate to it

Step 2: Download test graph and ch (area of Stuttgart, Germany)
	- wget https://fmi.uni-stuttgart.de/files/alg/data/graphs/stgtregbz.fmi.bz2
	- wget https://fmi.uni-stuttgart.de/files/alg/data/graphs/stgtregbz_ch.fmi.bz2

Step 3: Extract files
	- bzip2 -d stgtregbz.fmi.bz2
	- bzip2 -d stgtregbz_ch.fmi.bz2

Step 4: Compile code
	- g++ -O3 ch-light.cpp -o ch-light

Step 5: Run code
	- ./ch-light graphs/stgtregbz_ch.fmi graphs/stgtregbz.fmi

