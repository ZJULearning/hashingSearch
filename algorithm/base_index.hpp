#ifndef EFANNA_BASE_INDEX_H_
#define EFANNA_BASE_INDEX_H_

#include "general/params.hpp"
#include "general/distance.hpp"
#include "general/matrix.hpp"
#include <boost/dynamic_bitset.hpp>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <unordered_set>
#include <mutex>
#include "boost/smart_ptr/detail/spinlock.hpp"
#include <memory>



//#define BATCH_SIZE 200

namespace efanna{

typedef boost::detail::spinlock Lock;
typedef std::lock_guard<Lock> LockGuard;
struct Point {
	unsigned id;
	float dist;
	bool flag;
	Point () {}
	Point (unsigned i, float d, bool f = true): id(i), dist(d), flag(f) {
	}
	bool operator < (const Point &n) const{
		return this->dist < n.dist;
	}
};

typedef std::vector<Point>  Points;
static inline unsigned InsertIntoKnn (Point *addr, unsigned K, Point nn) {
	// find the location to insert
	unsigned j;
	unsigned i = K;
	while (i > 0) {
		j = i - 1;
		if (addr[j].dist <= nn.dist) break;
		i = j;
	}
	// check for equal ID
	unsigned l = i;
	while (l > 0) {
		j = l - 1;
		if (addr[j].dist < nn.dist) break;
		if (addr[j].id == nn.id) return K + 1;
		l = j;
	}
	// i <= K-1
	j = K;
	while (j > i) {
		addr[j] = addr[j-1];
		--j;
	}
	addr[i] = nn;
	return i;
}
struct Neighbor {
	std::shared_ptr<Lock> lock;
	float radius;
	float radiusM;
	Points pool;
	unsigned L;
	unsigned Range;
	bool found;
	std::vector<unsigned> nn_old;
	std::vector<unsigned> nn_new;
	std::vector<unsigned> rnn_old;
	std::vector<unsigned> rnn_new;

	Neighbor() : lock(std::make_shared<Lock>())
	{
	}

	unsigned insert (unsigned id, float dist) {
		if (dist > radius) return pool.size();
		LockGuard guard(*lock);
		unsigned l = InsertIntoKnn(&pool[0], L, Point(id, dist, true));
		if (l <= L) {
			if (L + 1 < pool.size()) {
				++L;
			}
			else {
				radius = pool[L-1].dist;
			}
		}
		return l;
	}

	template <typename C>
	void join (C callback) const {
		for (unsigned const i: nn_new) {
			for (unsigned const j: nn_new) {
				if (i < j) {
					callback(i, j);
				}
			}
			for (unsigned j: nn_old) {
				callback(i, j);
			}
		}
	}

};

template <typename DataType>
class InitIndex{
public:
	InitIndex(const Matrix<DataType>& features, const Distance<DataType>* d, const IndexParams& params):
		features_(features),
		distance_(d),
		params_(params)
{
}
	virtual ~InitIndex() {};
	virtual void buildTrees(){}

	virtual void buildIndex()
	{
		buildIndexImpl();
	}
	virtual void buildIndexImpl() = 0;
	virtual void loadGraph(char* filename) = 0;
	virtual void outputVisitBucketNum() = 0;

	void saveResults(char* filename){
		std::ofstream out(filename,std::ios::binary);
		std::vector<std::vector<int>>::iterator i;
		//std::cout<<nn_results.size()<<std::endl;
		for(i = nn_results.begin(); i!= nn_results.end(); i++){
			std::vector<int>::iterator j;
			int dim = i->size();
			//std::cout<<dim<<std::endl;
			out.write((char*)&dim, sizeof(int));
			for(j = i->begin(); j != i->end(); j++){
				int id = *j;
				out.write((char*)&id, sizeof(int));
			}
		}
		out.close();
	}
	SearchParams SP;
	void setSearchParams(int epochs, int init_num, int extend_to, int search_method){
		SP.search_epoches = epochs;
		SP.search_init_num = init_num;
		if(extend_to>init_num) SP.extend_to = init_num;
		else  SP.extend_to = extend_to;
		SP.search_method = search_method;
	}


	void nnExpansion_kgraph(size_t K, const DataType* qNow, std::vector<std::pair<float,size_t>>& pool, std::vector<int>& res, bool bSorted){
		unsigned int base_n = features_.get_rows();
		boost::dynamic_bitset<> tbflag(base_n, false);

		std::vector<Point> knn(K + SP.extend_to +1);
		std::vector<Point> results;

		for (unsigned iter = 0; iter < (unsigned)SP.search_epoches; iter++) {

			unsigned L = 0;
			for(unsigned j=0; j < (unsigned)SP.extend_to ; j++){
				if(!tbflag.test(pool[iter*SP.extend_to+j].second)){
					tbflag.set(pool[iter*SP.extend_to+j].second);
					knn[L].id = pool[iter*SP.extend_to+j].second;
					knn[L].dist = pool[iter*SP.extend_to+j].first;
					knn[L].flag = true;
					L++;
				}
			}
			if(~bSorted){
				std::sort(knn.begin(), knn.begin() + L);
			}

			unsigned int k =  0;
			while (k < L) {
				unsigned int nk = L;
				if (knn[k].flag) {
					knn[k].flag = false;
					unsigned n = knn[k].id;

					for(unsigned neighbor=0; neighbor < gs[n].size(); neighbor++){
						unsigned id = gs[n][neighbor];

						if(tbflag.test(id))continue;
						tbflag.set(id);

						float dist = distance_->compare(qNow, features_.get_row(id), features_.get_cols());
						Point nn(id, dist);
						unsigned int r = InsertIntoKnn(&knn[0], L, nn);

						//if ( (r <= L) && (L + 1 < knn.size())) ++L;
						if ( L + 1 < knn.size()) ++L;
						if (r < nk) nk = r;
					}
				}
				if (nk <= k) k = nk;
				else ++k;
			}


			if (L > K) L = K;
			if (results.empty()) {
				results.reserve(K + 1);
				results.resize(L + 1);
				std::copy(knn.begin(), knn.begin() + L, results.begin());
			}
			else {
				for (unsigned int l = 0; l < L; ++l) {
					unsigned r = InsertIntoKnn(&results[0], results.size() - 1, knn[l]);
					if (r < results.size() /* inserted */ && results.size() < (K + 1)) {
						results.resize(results.size() + 1);
					}
				}
			}
		}

		for(size_t i = 0; i < K && i < results.size();i++)
			res.push_back(results[i].id);
	}



	void nnExpansion(size_t K, const DataType* qNow, std::vector<int>& pool, boost::dynamic_bitset<>& tbflag){

		int resultSize = SP.extend_to;
		if (K > (unsigned)SP.extend_to)
			resultSize = K;

		unsigned int base_n = features_.get_rows();
		boost::dynamic_bitset<> newflag(base_n, true);
		newflag.set();

		std::vector<std::pair<float,size_t>> result;
		for(unsigned int i=0; i<pool.size();i++){
			result.push_back(std::make_pair(distance_->compare(qNow, features_.get_row(pool[i]), features_.get_cols()),pool[i]));
		}
		std::partial_sort(result.begin(), result.begin() + resultSize, result.end());
		result.resize(resultSize);
		pool.clear();
		for(int j = 0; j < resultSize; j++)
			pool.push_back(result[j].second);

		int iter=0;
		std::vector<int> ids;
		while(iter++ < SP.search_epoches){
			//the heap is max heap
			ids.clear();
			for(unsigned j = 0; j < SP.extend_to ; j++){
				if(newflag.test( pool[j] )){
					newflag.reset(pool[j]);

					for(unsigned neighbor=0; neighbor < gs[pool[j]].size(); neighbor++){
						unsigned id = gs[pool[j]][neighbor];

						if(tbflag.test(id))continue;
						else tbflag.set(id);

						ids.push_back(id);
					}
				}
			}

			for(size_t j = 0; j < ids.size(); j++){
				result.push_back(std::make_pair(distance_->compare(qNow, features_.get_row(ids[j]), features_.get_cols()),ids[j]));
			}
			std::partial_sort(result.begin(), result.begin() + resultSize, result.end());
			result.resize(resultSize);
			pool.clear();
			for(int j = 0; j < resultSize; j++)
				pool.push_back(result[j].second);
		}

		if(K<(unsigned)SP.extend_to)
			pool.resize(K);

	}

	virtual void knnSearch(int K, const Matrix<DataType>& query){
		getNeighbors(K,query);
	}

	virtual void getNeighbors(size_t K, const Matrix<DataType>& query) = 0;

	virtual void locatPoints(int nT, const Matrix<DataType>& query){
		locateNeighbors(nT,query);
	}

	virtual void locateNeighbors(size_t nT, const Matrix<DataType>& query) = 0;


protected:
	const Matrix<DataType> features_;
	const Distance<DataType>* distance_;
	const IndexParams params_;
	std::vector<std::vector<unsigned> > gs;
	std::vector<std::vector<int> > nn_results;
};
#define USING_BASECLASS_SYMBOLS \
		using InitIndex<DataType>::distance_;\
		using InitIndex<DataType>::params_;\
		using InitIndex<DataType>::features_;\
		using InitIndex<DataType>::buildIndex;\
		using InitIndex<DataType>::nn_results;\
		using InitIndex<DataType>::saveResults;\
		using InitIndex<DataType>::SP;\
		using InitIndex<DataType>::nnExpansion;\
		using InitIndex<DataType>::nnExpansion_kgraph;\
		using InitIndex<DataType>::gs;\

}
#endif
