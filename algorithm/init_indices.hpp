#ifndef EFANNA_INIT_INDICES_H_
#define EFANNA_INIT_INDICES_H_
#include "base_index.hpp"
#include "hashing_index.hpp"
#include "hamming_index.hpp"
#include "ikhamming_index.hpp"

namespace efanna{


  template <template<typename> class Index, typename DataType>
  inline InitIndex<DataType>* create_index_(efanna::Matrix<DataType> data, const efanna::IndexParams& params,  const Distance<DataType>* d)
  {
      return new Index<DataType>(data, d, params);
  }

  template <typename DataType>
  inline InitIndex<DataType>*
  create_index_by_type(const init_algorithm index_type,
    		const Matrix<DataType>& dataset, const IndexParams& params,  const Distance<DataType>* d)
  {
	  InitIndex<DataType>* initIndex = NULL;
	  switch(index_type){
	  case HAMMING:
		  initIndex = create_index_<HAMMINGIndex, DataType>(dataset, params, d);
		  break;
	  case HASHING:
		  initIndex = create_index_<HASHINGIndex, DataType>(dataset, params, d);
		  break;
      case IKHAMMING:
    	  initIndex = create_index_<IKHAMMINGIndex, DataType>(dataset, params, d);
    	  break;
	  }
	  return initIndex;
  }
}
#endif
