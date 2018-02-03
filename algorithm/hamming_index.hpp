#ifndef EFANNA_HAMMING_INDEX_H_
#define EFANNA_HAMMING_INDEX_H_

#include "algorithm/base_index.hpp"
#include <boost/dynamic_bitset.hpp>
#include <time.h>
//for Debug
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <map>
#include <sstream>
#include <set>



namespace efanna{


inline unsigned parallel_popcnt32(unsigned x) {
	x = (x & 0x55555555) + ((x >> 1 ) & 0x55555555);
	x = (x & 0x33333333) + ((x >> 2 ) & 0x33333333);
	x = (x & 0x0f0f0f0f) + ((x >> 4 ) & 0x0f0f0f0f);
	x = (x & 0x00ff00ff) + ((x >> 8 ) & 0x00ff00ff);
	x = (x & 0x0000ffff) + ((x >> 16) & 0x0000ffff);
	return x;
}

struct HAMMINGIndexParams : public IndexParams
{
	HAMMINGIndexParams(int codelen, char*& BaseCodeFile, char*& QueryCodeFile)
	{
		init_index_type = HAMMING;
		ValueType len;
		len.int_val = codelen;
		extra_params.insert(std::make_pair("codelen",len));
		ValueType bcf;
		bcf.str_pt = BaseCodeFile;
		extra_params.insert(std::make_pair("bcfile",bcf));
		ValueType qcf;
		qcf.str_pt = QueryCodeFile;
		extra_params.insert(std::make_pair("qcfile",qcf));
	}
};

template <typename DataType>
class HAMMINGIndex : public InitIndex<DataType>
{
public:

	typedef InitIndex<DataType> BaseClass;

	typedef std::vector<unsigned int> Code;
	typedef std::vector<Code> Codes;


	HAMMINGIndex(const Matrix<DataType>& dataset, const Distance<DataType>* d, const IndexParams& params = HAMMINGIndexParams(0,NULL,NULL)) :
		BaseClass(dataset,d,params)
	{
		std::cout<<"HASHING initial" <<std::endl;

		ExtraParamsMap::const_iterator it = params_.extra_params.find("codelen");
		if(it != params_.extra_params.end()){
			codelength = (it->second).int_val;
		}
		else{
			std::cout << "error: no code length setting" << std::endl;
		}

		it = params_.extra_params.find("bcfile");
		if(it != params_.extra_params.end()){
			char* fpath = (it->second).str_pt;
			std::string str(fpath);
			std::cout << "Loading base code from " << str << std::endl;

			LoadCode(fpath, BaseCode);
		}
		else{
			std::cout << "error: no base code file" << std::endl;
		}

		it = params_.extra_params.find("qcfile");
		if(it != params_.extra_params.end()){
			char* fpath = (it->second).str_pt;
			std::string str(fpath);
			std::cout << "Loading query code from " << str << std::endl;

			LoadCode(fpath, QueryCode);

		}
		else{
			std::cout << "error: no query code file" << std::endl;
		}

		std::cout << "Search begin..." << std::endl;
	}


	void LoadCode(char* filename, Codes& base){

		unsigned tableNum = codelength / 32;
		unsigned lastLen = codelength % 32;

		std::stringstream ss;
		unsigned j;
		for(j = 0; j < tableNum; j++){

			ss << filename << "_" << j+1 ;
			std::string codeFile;
			ss >> codeFile;
			ss.clear();

			std::ifstream in(codeFile.c_str(), std::ios::binary);
			if(!in.is_open()){std::cout<<"open file " << filename <<" error"<< std::endl;return;}

			int codeNum;
			in.read((char*)&codeNum,4);
			if (codeNum != 1){
				std::cout<<"codefile  "<< j << " error!"<<std::endl;
			}

			unsigned len;
			in.read((char*)&len,4);
			unsigned num;
			in.read((char*)&num,4);

			Code b;
			for(unsigned i = 0; i < num; i++){
				unsigned int codetmp;
				in.read((char*)&codetmp,4);
				b.push_back(codetmp);
			}
			base.push_back(b);
			in.close();
		}

		if(lastLen > 0){
			int shift = 32 - lastLen;
			ss << filename << "_" << j+1 ;
			std::string codeFile;
			ss >> codeFile;
			ss.clear();

			std::ifstream in(codeFile.c_str(), std::ios::binary);
			if(!in.is_open()){std::cout<<"open file " << filename <<" error"<< std::endl;return;}

			int codeNum;
			in.read((char*)&codeNum,4);
			if (codeNum != 1){
				std::cout<<"codefile  "<< j << " error!"<<std::endl;
			}

			unsigned len;
			in.read((char*)&len,4);
			unsigned num;
			in.read((char*)&num,4);

			Code b;
			for(unsigned i = 0; i < num; i++){
				unsigned int codetmp;
				in.read((char*)&codetmp,4);
				codetmp = codetmp >> shift;
				b.push_back(codetmp);
			}
			base.push_back(b);
			in.close();
		}

	}


	void buildIndexImpl(){}


	void locateNeighbors(size_t nT, const Matrix<DataType>& query){

		size_t pNum = BaseCode[0].size();
		size_t nTable = BaseCode.size();

		//if (nT > nTable){std::cout << "index code length not enough!"; nT = nTable;}
		nT = nTable;

		for (size_t cur = 0; cur < query.get_rows(); cur++) {
			std::vector<int>hammingDistance(pNum);
			std::fill(hammingDistance.begin(), hammingDistance.end(), 0);

			for (size_t j = 0; j < nT; j++) {
				for (size_t i = 0; i < pNum; i++) {
					hammingDistance[i] += parallel_popcnt32(BaseCode[j][i] ^ QueryCode[j][cur]);
				}
			}

			// initialize index
			std::vector<size_t> hammingDistanceIndex(pNum);
			for (size_t i = 0; i < hammingDistanceIndex.size(); i++) {
				hammingDistanceIndex[i] = i;
			}
			// compare function
			auto compare_func = [&hammingDistance](const size_t a, const size_t b)->bool
					{
				return hammingDistance[a] < hammingDistance[b];
					};
			// sort the index by its value
			std::partial_sort(hammingDistanceIndex.begin(), hammingDistanceIndex.begin() + SP.search_init_num, hammingDistanceIndex.end(),  compare_func);

		}

	}


	void getNeighbors(size_t K, const Matrix<DataType>& query){

		size_t pNum = BaseCode[0].size();
		size_t nTable = BaseCode.size();

		nn_results.clear();
		for (size_t cur = 0; cur < query.get_rows(); cur++) {
			std::vector<int>hammingDistance(pNum);
			std::fill(hammingDistance.begin(), hammingDistance.end(), 0);

			for (size_t j = 0; j < nTable; j++) {
				for (size_t i = 0; i < pNum; i++) {
					hammingDistance[i] += parallel_popcnt32(BaseCode[j][i] ^ QueryCode[j][cur]);
				}
			}

			// initialize index
			std::vector<size_t> hammingDistanceIndex(pNum);
			for (size_t i = 0; i < hammingDistanceIndex.size(); i++) {
				hammingDistanceIndex[i] = i;
			}
			// compare function
			auto compare_func = [&hammingDistance](const size_t a, const size_t b)->bool
					{
				return hammingDistance[a] < hammingDistance[b];
					};
			// sort the index by its value
			std::partial_sort(hammingDistanceIndex.begin(), hammingDistanceIndex.begin() + SP.search_init_num, hammingDistanceIndex.end(),  compare_func);

			std::vector<int> res;
			if ((unsigned) SP.search_init_num <= K){
				for (unsigned int j = 0; j < (unsigned) SP.search_init_num; j++) res.push_back(hammingDistanceIndex[j]);
			}else{
				std::vector<std::pair<float, size_t>> Distance;
				for (size_t i = 0; i < (unsigned) SP.search_init_num; i++) {
					Distance.push_back(std::make_pair(distance_->compare(query.get_row(cur), features_.get_row(hammingDistanceIndex[i]), features_.get_cols()), hammingDistanceIndex[i]));
				}
				std::partial_sort(Distance.begin(), Distance.begin() + K, Distance.end());

				for (unsigned int j = 0; j < K; j++) res.push_back(Distance[j].second);
				nn_results.push_back(res);
			}
		}
	}






	void outputVisitBucketNum(){
	}


	void loadGraph(char* filename){
		std::ifstream in(filename,std::ios::binary);
		in.seekg(0,std::ios::end);
		std::ios::pos_type ss = in.tellg();
		size_t fsize = (size_t)ss;
		int dim;
		in.seekg(0,std::ios::beg);
		in.read((char*)&dim, sizeof(int));
		size_t num = fsize / (dim+1) / 4;
		std::cout<<"load g "<<num<<" "<<dim<< std::endl;
		in.seekg(0,std::ios::beg);
		gs.clear();
		for(size_t i = 0; i < num; i++){
			std::vector<unsigned> heap;
			in.read((char*)&dim, sizeof(int));
			for(int j =0; j < dim; j++){
				unsigned id;
				in.read((char*)&id, sizeof(int));
				heap.push_back(id);
			}
			gs.push_back(heap);
		}
		in.close();
	}


protected:
	USING_BASECLASS_SYMBOLS

	int codelength;
	Codes BaseCode;
	Codes QueryCode;

};
}
#endif
