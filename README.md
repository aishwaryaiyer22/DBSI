# DBSI
	W4112 Database Systems Implementation
		     Project 2


Ranjith Kumar Shanmuga Visvanathan	- rs3579
Aishwarya Iyer 				- ai2336


Part 1
======

	- Pointer free B-Tree implementation in C
	- Number of keys and number of probes obtained as arguments
	- Also the fanouts in each level is obtained as arguments
		Eg. build 10 12 3 3 3

	
	Implementation :
	----------------

		- First the array of array pointer, which maintains the pointer to the array of the each level of the tree. The first element points to the level 0 array of the tree, and the next to level 1 and so on.
		
		- Memory allocation for each level is done based on the fanout and size of each level. Each level is allocated the max size.

		- The size for each level is calcuated based on the size of the previous level, from which the number of nodes/blocks for this level is figured. This times the size of each block gives  the total size of this level.

		- Also posix_memalign with alignment of 16 is used to allocate the memory of each level. So they all start at addresses multiple of 16.

		- Now with the total tree size known, we can say if the keys are more than tree size.

		 Insertion:
		 -----------
		
		- Bulk Insertion method is followed, where we start inserting from the bottom most node, and whenever each node at a level is full, we move to its parent. If it happes there we move to its parent and so on.

		- this is accomplished by keeping reference variables of current position for each level.

		 Key Generation & Sorting:
		 --------------------------

		- The reference code "p2random v2.c" is converted to a headerfile "p2randomv2.h" and used
		- The functions generate_sorted_unique is used to get a random set of the number of keys requested and also sort them and returnes the pointer to the array. WHich is inserted into the B Tree.

		 Probe Generation:
		 ------------------
		
		- The probe generation also uses the generate function from the given reference code

		
		
		Range Identifier:
		----------------

		- We start from the root and use binary search, a recursive binary search is implemented, to figure out which range identifier its falling in each level. Its then added to the array, and we move to the next level below by identifying the correct block.
		- Identifying the block is important as its needed in giving the start and end positions for the binary search for that level.
		- The start for the next level can be calucated from the current level, current block and the range identifier of the current block
		
		- Start is basically the position of the node block we are going to binary search next. This is calculated from the previous row block postion (like the oth/1st/3rd node block etc) and the fanout of the previous level
		
		- Once start is figured end is start + the size of block - 1

		- Again BS can be called to get the range identifier for this level and then get the range identifier for this level, and so on until the last level.

		- At each level we know, whats the starting position, and also the range identification for that node, so with this, we can canculate the total number of nodes before at each level. This cumulative sum gives us the Range identification Thus we can arrive at the overal Range Identifier for the given probe.

		- Hence we have the Range Identification for the given probe, also the identification at each node level on the path from root to leaf too..
