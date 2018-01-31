Search with A Hash Index 
============
We provide here the codes for search with a hash index (and a kNN graph). The codes are used in our paper [A Revisit of Hashing Algorithms for Approximate Nearest Neighbor Search](http://arxiv.org/abs/1612.07545).    

You need a hash index for search. Please use the [matlab codes](https://github.com/dengcai78/MatlabFunc/tree/master/ANNS/Hashing) of various hashing algorithms to generate the hash index.

Benchmark data set
-------
* [SIFT1M and GIST1M](http://corpus-texmex.irisa.fr/)

The performance was tested without parallelism.   

ANN search using hash index (hamming ranking approach)
------

![SIFT100nn](figures/SIFT.png)     
![GIST100nn](figures/GIST.png)    

The number in parenthesis after the algorithm name is the optimal code length of this hashing algorithm on this dataset.


How To Complie    
-------
Go to the root directory of EFANNA and make.    

	cd hashingSearch/
	make

How To Use    
------
You need a hash index for search. Please use the [matlab codes](https://github.com/dengcai78/MatlabFunc/tree/master/ANNS/Hashing) of various hashing algorithms to generate the hash index.

* ANN search with a hash index using hamming ranking approach

		cd samples/
		./hamming_search  data_file BaseCodeFile query_file QueryCodeFile result_file codelen initsz querNN

  Meaning of the parameters(from left to right):   

	data_file     -- database points  
	BaseCodeFile  -- e.g. LSHtableSIFT32b. (The actual binary code file is LSHtableSIFT32b_1, LSHtableSIFT32b_2, ..., the program will automatically loaded the required number of bits specified by the codelen parameter)  
	query_file    -- sift query points  
	QueryCodeFile -- e.g. LSHquerySIFT32b. (The actual binary code file is LSHquerySIFT32b_1, LSHquerySIFT32b_2, ..., the program will automatically loaded the required number of bits specified by the codelen parameter)  
	result_file   -- path to save ANN search results of given query   
	codelen       -- code length of the binary codes. (A long code will be split as multiple 32-bits tables for convenience)   
	initsz        -- initial pool size, the parameter L in the paper.
	querNN        -- required number of returned neighbors (i.e. k of k-NN)   


* ANN search with a hash index using hash bucket search approach

		cd samples/
		./hashing_search sift_base.fvecs 16 LSHtableSIFT32b sift_query.fvecs LSHquerySIFT32b sift.results 16 32 8 0 200 100

  Meaning of the parameters(from left to right):   

	sift_base.fvecs -- database points  
	LSHtableSIFT32b -- the binary code file of the database (the actual file is LSHtableSIFT32b_1, LSHtableSIFT32b_2, ...)  
	sift_query -- sift query points  
	LSHquerySIFT32b -- the binary code file of the query (the actual file is LSHquerySIFT32b_1, LSHquerySIFT32b_2, ...)  
	sift.results -- path to save ANN search results of given query   
	16 -- number of tables to use (you should provide 16 tables)   
	32 -- code length of the hash table   
	8  -- maximum radius (the algorithm will stop locating points beyond this radius number, using random points from the database to fill the initial pool)   
	0  -- code length shift (the algorithm will actually use codelength - codelengthshift bits code, in this case, 32-0=32)   
    200 -- initial pool size factor    
	100 -- required number of returned neighbors (i.e. k of k-NN)   


Output and Input format
------
Same as that of [EFANNA](https://github.com/fc731097343/efanna)



