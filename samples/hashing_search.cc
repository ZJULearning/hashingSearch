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

  int align_dim = (dim + 7)/8*8;
  cout<<"data aligned dimension: "<<align_dim<<endl;

  data = new float[num*align_dim]();

  in.seekg(0,ios::beg);
  for(size_t i = 0; i < num; i++){
    in.seekg(4,ios::cur);
    in.read((char*)(data+i*align_dim),dim*4);
  }
  in.close();

  dim = align_dim;
}
int main(int argc, char** argv){
	if(argc!=10 && argc!=11){cout<< argv[0] << " data_file BaseCodeFile query_file QueryCodeFile result_file tablelen codelen initsz querNN index_method(optional)" <<endl; exit(-1);}

  float* data_load = NULL;
  float* query_load = NULL;
  size_t points_num;
  int dim;
  load_data(argv[1], data_load, points_num,dim);
  size_t q_num;
  int qdim;
  load_data(argv[3], query_load, q_num,qdim);
  Matrix<float> dataset(points_num,dim,data_load);
  Matrix<float> query(q_num,qdim,query_load);

  int tablelen = atoi(argv[6]);
  int codelen = atoi(argv[7]);
  int index_method = argc == 11 ? atoi(argv[10]) : 1;
  FIndex<float> index(dataset, new L2DistanceAVX<float>(), efanna::HASHINGIndexParams(codelen, tablelen,argv[2],argv[4],index_method));

  auto s = std::chrono::high_resolution_clock::now();
  index.buildIndex();
  auto e = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = e-s;
  std::cout << "Index building time: " << diff.count() << "\n";

  int poolsz = atoi(argv[8]);
  index.setSearchParams(codelen, poolsz, poolsz);


  s = std::chrono::high_resolution_clock::now();
  index.knnSearch(atoi(argv[9])/*query nn*/,query);
  e = std::chrono::high_resolution_clock::now();
  diff = e-s;
  std::cout << "query searching time: " << diff.count() << "\n";


  index.saveResults(argv[5]);

  index.outputVisitBucketNum();

  return 0;
}
