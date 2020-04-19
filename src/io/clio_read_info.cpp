////////////////////////////////////////////////////////////////////////
// IO: basics
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
////////////////////////////////////////////////////////////////////////
void clio::read_info(std::string key, std::string& val) const
{
  /* Read a bool val labeled by key from inp. */

  std::ifstream _ifs( inp.c_str() );
  _ifs.seekg( 0, std::ios::beg );

  std::string line;
  std::stringstream ioss;

  key.append(" ");
  long ind_key;

  //  while ( getline(_ifs, line) && line.find(key,0) == std::string::npos ) {}
  while ( getline(_ifs, line) && line.substr(0, key.length()) != key ) {}

  if(! _ifs.eof()) {

    ind_key = line.find( "=" ) + 1;

    ioss.str("");
    ioss << line.substr( ind_key );
    ioss >> val;

    std::cout << "# input-ess " << key.c_str() << " = " << val << std::endl;

  } else {

    std::cout << "# " << "Input error:" << key << " not found." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clio::read_info(std::string key, bool& val) const
{
  /* Read a bool val labeled by key from inp. */

  std::ifstream _ifs( inp.c_str() );
  _ifs.seekg( 0, std::ios::beg );

  std::string line;
  std::string val_string;
  std::stringstream ioss;

  key.append(" ");
  long ind_key;

  //  while ( getline(_ifs, line) && line.find(key,0) == std::string::npos ) {}
  while ( getline(_ifs, line) && line.substr(0, key.length()) != key ) {}

  if(! _ifs.eof()) {

    ind_key = line.find( "=" ) + 1;

    ioss.str("");
    ioss << line.substr( ind_key );
    ioss >> val_string;

    std::cout << "# input-ess " << key.c_str() << " = " << val_string << std::endl;

  } else {

    std::cout << "# " << "Input error:" << key << " not found." << std::endl;
    abort();
  }

  if (val_string.compare("true") == 0) {
    val = true;
  } else if (val_string.compare("false") == 0) {
    val = false;
  } else {
    std::cout << "Input error: set true or false for " << key << "." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clio::read_info(std::string key, long& val) const
{
  /* Read an integer val labeled by key from inp. */

  std::ifstream _ifs( inp.c_str() );
  _ifs.seekg( 0, std::ios::beg );

  std::string line;
  std::stringstream ioss;

  key.append(" ");
  long ind_key;

  //  while ( getline(_ifs, line) && line.find(key,0) == std::string::npos ) {}
  while ( getline(_ifs, line) && line.substr(0, key.length()) != key ) {}

  if(! _ifs.eof()) {

    ind_key = line.find( "=" ) + 1;

    ioss.str("");
    ioss << line.substr( ind_key );
    ioss >> val;

    std::cout << "# input-ess " << key.c_str() << " = " << val << std::endl;

  } else {

    std::cout << "# " << "Input error:" << key << " not found." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clio::read_info(std::string key, double& val) const
{
  /* Read an double val labeled by key from inp. */

  std::ifstream _ifs( inp.c_str() );
  _ifs.seekg( 0, std::ios::beg );

  std::string line;
  std::stringstream ioss;

  key.append(" ");
  long ind_key;

  //  while ( getline(_ifs, line) && line.find(key,0) == std::string::npos ) {}
  while ( getline(_ifs, line) && line.substr(0, key.length()) != key ) {}

  if(! _ifs.eof()) {

    ind_key = line.find( "=" ) + 1;

    ioss.str("");
    ioss << line.substr( ind_key );
    ioss >> val;

    std::cout << "# input-ess " << key.c_str() << " = " << val << std::endl;

  } else {

    std::cout << "Input error:" << key << " not found." << std::endl;
    abort();

  }

}
////////////////////////////////////////////////////////////////////////
void clio::read_info(std::string key, dcomplex& val) const
{
  /* Read an complex val labeled by key from inp. */

  std::ifstream _ifs( inp.c_str() );
  _ifs.seekg( 0, std::ios::beg );

  std::string line;
  std::stringstream ioss;

  key.append(" ");
  long ind_key;

  //  while ( getline(_ifs, line) && line.find(key,0) == std::string::npos ) {}
  while ( getline(_ifs, line) && line.substr(0, key.length()) != key ) {}

  if(! _ifs.eof()) {

    ind_key = line.find( "=" ) + 1;

    ioss.str("");
    ioss << line.substr( ind_key );
    ioss >> val;

    std::cout << "# input-ess " << key.c_str() << " = " << val << std::endl;

  } else {

    std::cout << "Input error:" << key << " not found." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clio::read_info(std::string key, std::vector<long>& vals) const
{
  /* Read integer vals labeled by key from inp. */

  std::ifstream _ifs( inp.c_str() );
  _ifs.seekg( 0, std::ios::beg );

  std::string line, trash;
  std::stringstream ioss;

  key.append(" ");
  long ind_key;

  //DEBUG
  //  std::cout << "size = " << vals.size() << std::endl;
  //DEBUG

  //  while ( getline(_ifs, line) && line.find(key,0) == std::string::npos ) {}
  while ( getline(_ifs, line) && line.substr(0, key.length()) != key ) {}

  if(! _ifs.eof()) {

    line = line.substr(line.find("=") + 1);

    ioss.str("");
    ioss << line;
    for (long id = 0; id < vals.size(); id ++) {
      ioss >> vals[id];
      ioss >> trash;
    }

//    line = line.substr(line.find("=") + 1);
//    for (long id = 0; id < vals.size(); id ++) {
//      ioss.str("");
//      ioss << line.substr(line.find(","));
//      ioss >> vals[id];
//      line = line.substr(line.find(",")+1);
//    }
//
    std::cout << "# input-ess " << key.c_str() << " = ";    
    for (long id = 0; id < vals.size() - 1; id ++) std::cout << vals[id] << ", ";
    std::cout << vals[vals.size()-1] << std::endl;

  } else {

    std::cout << "# " << "Input error:" << key << " not found." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
