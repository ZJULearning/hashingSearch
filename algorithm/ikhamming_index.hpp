#ifndef EFANNA_IKHAMMING_INDEX_H_
#define EFANNA_IKHAMMING_INDEX_H_

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
#include <numeric>
#include <utility>



namespace efanna{



struct IKHAMMINGIndexParams : public IndexParams
{
	IKHAMMINGIndexParams(int codelen, char*& BaseCodeFile, char*& QueryCodeFile, char*& InvertedKmeansFile)
	{
		init_index_type = IKHAMMING;
		ValueType len;
		len.int_val = codelen;
		extra_params.insert(std::make_pair("codelen",len));
		ValueType bcf;
		bcf.str_pt = BaseCodeFile;
		extra_params.insert(std::make_pair("bcfile",bcf));
		ValueType qcf;
		qcf.str_pt = QueryCodeFile;
		extra_params.insert(std::make_pair("qcfile",qcf));
		ValueType ikf;
		ikf.str_pt = InvertedKmeansFile;
		extra_params.insert(std::make_pair("ikfile",ikf));
	}
};

template <typename DataType>
class IKHAMMINGIndex : public InitIndex<DataType>
{
public:

	typedef InitIndex<DataType> BaseClass;

	typedef std::vector<unsigned int> Code;
	typedef std::vector<Code> Codes;
	typedef std::vector<Codes> BucketCodes;


	IKHAMMINGIndex(const Matrix<DataType>& dataset, const Distance<DataType>* d, const IndexParams& params = IKHAMMINGIndexParams(0,NULL,NULL,NULL)) :
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

		it = params_.extra_params.find("ikfile");
		if(it != params_.extra_params.end()){
			char* fpath = (it->second).str_pt;
			std::string str(fpath);
			std::cout << "Loading kmenas partitions from " << str << std::endl;

			LoadPartition(fpath);

		}
		else{
			std::cout << "error: no kmeans file" << std::endl;
		}
	}


	void LoadCode(char* filename, Codes& base){

		unsigned tableNum = codelength / 32;
		unsigned lastLen = codelength % 32;
		if(lastLen > 0){
			std::cout <<"codelength not supported!" <<std::endl;
			return;
		}

		std::stringstream ss;
		for(unsigned j = 0; j < tableNum; j++){

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

	}

	void LoadPartition(char* filename){

		std::ifstream in(filename, std::ios::binary);
		if(!in.is_open()){std::cout<<"open file " << filename <<" error"<< std::endl;return;}


		unsigned nCluster;
		in.read((char*)&nCluster,4);

		unsigned Dim;
		in.read((char*)&Dim,4);

		centers.resize(nCluster);
		hashTable.resize(nCluster);

		unsigned tableNum = BaseCode.size();

		for(unsigned i = 0; i < nCluster; i++){

			Codes tmp(tableNum);
			BaseBucket.push_back(tmp);

			unsigned nSmp;
			in.read((char*)&nSmp,4);

			for(unsigned j = 0; j < Dim; j++){
				float tmp;
				in.read((char*)&tmp,4);
				centers[i].push_back(tmp);
			}

			for(unsigned j = 0; j < nSmp; j++){
				unsigned tmp;
				in.read((char*)&tmp,4);
				hashTable[i].push_back(tmp);

				for(unsigned k=0; k<tableNum; k++){
					BaseBucket[i][k].push_back(BaseCode[k][tmp]);
				}

			}


		}

		in.close();


	}


	void buildIndexImpl(){
		nGroup = SP.search_epoches;

		if(nGroup > centers.size()){
			nGroup = centers.size();
		}

		std::cout << "Search begin..." << std::endl;
	}



	void getNeighbors(size_t K, const Matrix<DataType>& query){


		size_t nTable = BaseCode.size();

		unsigned nCluster = centers.size();


		std::vector<std::vector<int> >hammingDistance(nGroup);
		std::vector<size_t> hammingDistanceIndex;;
		std::vector<std::pair<int, unsigned> >hammingDistanceSort;

		nn_results.clear();
		for (size_t cur = 0; cur < query.get_rows(); cur++) {


			std::vector<std::pair<float, unsigned> > CenterDistance;
			for(unsigned j = 0; j < nCluster; j++){
				CenterDistance.push_back(std::make_pair(distance_->compare(query.get_row(cur), centers[j].data(), query.get_cols()), j));
			}
			std::partial_sort(CenterDistance.begin(), CenterDistance.begin() + nGroup, CenterDistance.end());

			size_t pNum = 0;
			std::vector<unsigned> compareIndex;

			for(unsigned j = 0; j < nGroup; j++){
				compareIndex.push_back(CenterDistance[j].second);
				pNum += hashTable[CenterDistance[j].second].size();
			}
			if ((unsigned)SP.search_init_num > pNum){
				SP.search_init_num = pNum;
			}


			// resize hammingDistance
			for (size_t i = 0; i < compareIndex.size(); i++) {
				hammingDistance[i].resize(hashTable[compareIndex[i]].size());
				std::fill(hammingDistance[i].begin(), hammingDistance[i].end(), 0);
			}

			unsigned IndexNum = compareIndex.size();
			for (unsigned j = 0; j < nTable; j++) {
				for (unsigned idx = 0; idx < IndexNum; idx++) {
					unsigned tNum = BaseBucket[compareIndex[idx]][j].size();
					Code &code = BaseBucket[compareIndex[idx]][j];
					std::vector<int> &dst = hammingDistance[idx];
					for (unsigned i = 0; i < tNum; i++) {
						dst[i] += parallel_popcnt32(code[i] ^ QueryCode[j][cur]);
					}
				}
			}



			hammingDistanceSort.resize(0);
			for (unsigned idx = 0; idx < IndexNum; idx++) {
				size_t sort_num = hammingDistance[idx].size();
				for (size_t i = 0; i < sort_num; i++) {
					hammingDistanceSort.push_back(std::pair<int, unsigned>(hammingDistance[idx][i], hashTable[compareIndex[idx]][i]));
				}
			}



			// sort the index by its value
			size_t sort_num = std::min((size_t)SP.search_init_num, hammingDistanceSort.size());
			std::partial_sort(hammingDistanceSort.begin(), hammingDistanceSort.begin() + sort_num, hammingDistanceSort.end());

			std::vector<int> res;
			if ((unsigned) SP.search_init_num <= K){
				for (unsigned int j = 0; j < (unsigned) sort_num; j++) res.push_back(hammingDistanceSort[j].second);
			}else{
				std::vector<std::pair<float, unsigned> > Distance;
				for (size_t i = 0; i < (unsigned) sort_num; i++) {
					Distance.push_back(std::make_pair(distance_->compare(query.get_row(cur), features_.get_row(hammingDistanceSort[i].second), features_.get_cols()), hammingDistanceSort[i].second));
				}
				std::partial_sort(Distance.begin(), Distance.begin() + K, Distance.end());

				for (unsigned int j = 0; j < K; j++) res.push_back(Distance[j].second);
			}
			nn_results.push_back(res);
		}
	}


	void examineCodeImpl(){}
	void locateNeighbors(const Matrix<DataType>& query){}
	void outputVisitBucketNum(){}
	void loadGraph(char* filename){}


protected:
	USING_BASECLASS_SYMBOLS

	unsigned codelength;
	unsigned nGroup;

	Codes BaseCode;
	Codes QueryCode;

	std::vector<std::vector<float>> centers;
	std::vector<std::vector<unsigned>> hashTable;
	BucketCodes BaseBucket;

	// for statistic info
	std::vector<unsigned> HammingDistNum;

};
}
#endif
