#include <efanna.hpp>
#include <iostream>
#include <fstream>
#include <chrono>

using namespace efanna;
using namespace std;
void load_data(char* filename, float*& data, size_t& num,int& dim){// load data with sift10K pattern
  ifstream in(filename, ios::binary);
  if(!in.is_open()){cout<<"open file error"<<endl;exit(-1);}
  in.read((char*)&dim,4);
  cout<<"data dimension: "<<dim<<endl;
  in.seekg(0,ios::end);
  ios::pos_type ss = in.tellg();
  size_t fsize = (size_t)ss;
  num = fsize / (dim+1) / 4;
  data = new float[num*dim];

  in.seekg(0,ios::beg);
  for(size_t i = 0; i < num; i++){
    in.seekg(4,ios::cur);
    in.read((char*)(data+i*dim),dim*4);
  }
  in.close();
}
int main(int argc, char** argv){
	if(argc!=13 && argc!=16 && argc!=17){cout<< argv[0] << " data_file UpperBits BaseCodeFile query_file QueryCodeFile result_file tableNum codelen radius lenshift initsz  querNN graph_index(optional) epoch(optional) extend(optional) search_method(optional)" <<endl; exit(-1);}

  float* data_load = NULL;
  float* query_load = NULL;
  size_t points_num;
  int dim;
  load_data(argv[1], data_load, points_num,dim);
  size_t q_num;
  int qdim;
  load_data(argv[4], query_load, q_num,qdim);
  Matrix<float> dataset(points_num,dim,data_load);
  Matrix<float> query(q_num,qdim,query_load);

  int UpperBits = atoi(argv[2]);
  int tableNum = atoi(argv[7]);
  int codelen = atoi(argv[8]);
  int radius = atoi(argv[9]);
  int lenshift = atoi(argv[10]);
  FIndex<float> index(dataset, new L2DistanceAVX<float>(), efanna::HASHINGIndexParams(codelen, tableNum,UpperBits,radius,argv[3],argv[5],lenshift));

  auto s = std::chrono::high_resolution_clock::now();
  index.buildIndex();
  auto e = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = e-s;
  std::cout << "Index building time: " << diff.count() << "\n";

  int poolsz = atoi(argv[11]);
  if (argc == 13){
	  index.setSearchParams(tableNum, poolsz, poolsz);
  }else{
	  index.loadGraph(argv[13]);
	  int search_epoc = atoi(argv[14]);
	  int search_extend = atoi(argv[15]);
	  int search_method = argc == 17 ? atoi(argv[16]) : 0;
	  index.setSearchParams(search_epoc, poolsz, search_extend,search_method);
  }


  s = std::chrono::high_resolution_clock::now();
  index.knnSearch(atoi(argv[12])/*query nn*/,query);
  e = std::chrono::high_resolution_clock::now();
  diff = e-s;
  std::cout << "query searching time: " << diff.count() << "\n";


  index.saveResults(argv[6]);

  index.outputVisitBucketNum();

  return 0;
}
