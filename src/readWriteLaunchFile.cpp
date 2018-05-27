
// type support
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/complex.hpp>

// for doing the actual serialization
#include <cereal/archives/json.hpp>
#include <iostream>
#include <fstream>



struct data{
    int i1;
    int i2;
    
      template <class Archive>
  void serialize( Archive & ar )
  {
    ar( i1, i2 );
  }
};


void readLaunchFile(char *JsonFilename, data *pointerToData){
    data Data = *pointerToData;
    std::ifstream is(JsonFilename);
    cereal::JSONInputArchive ar(is);
    ar( cereal::make_nvp("var1", Data.i1) );
}

void writeLaunchFile(char *outFileName, data *pointerToData){
    data Data = *pointerToData;
    std:ofstream outFile(outFileName);
    cereal::JSONOutputArchive output(outFile);
    output( cereal::make_nvp("best data ever", Data) );
}
