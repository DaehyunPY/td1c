////////////////////////////////////////////////////////////////////////
// IO: basics
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
////////////////////////////////////////////////////////////////////////
void clio::read_info(std::string key, const std::string def, std::string& val) const
{
  /* Read a string val labeled by key from inp. */

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
  } else {
    val = def;
  }

  std::cout << "# input-opt " << key.c_str() << " = " << val << std::endl;
}
////////////////////////////////////////////////////////////////////////
void clio::read_info(std::string key, bool def, bool& val) const
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

    if (val_string.compare("true") == 0) {
      val = true;
    } else if (val_string.compare("false") == 0) {
      val = false;
    } else {
      std::cout << "Input error: set true or false for " << key << "." << std::endl;
      abort();
    }
  } else {
    val = def;
    if (val) {
      val_string = "true";
    } else {
      val_string = "false";
    }
  }

  std::cout << "# input-opt " << key.c_str() << " = " << val_string << std::endl;
}
////////////////////////////////////////////////////////////////////////
void clio::read_info(std::string key, long def, long& val) const
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
  } else {
    val = def;
  }

  std::cout << "# input-opt " << key.c_str() << " = " << val << std::endl;
}
////////////////////////////////////////////////////////////////////////
void clio::read_info(std::string key, double def, double& val) const
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
  } else {
    val = def;
  }
  
  std::cout << "# input-opt " << key.c_str() << " = " << val << std::endl;
}
////////////////////////////////////////////////////////////////////////
void clio::read_info(std::string key, dcomplex def, dcomplex& val) const
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
  } else {
    val = def;
  }

  std::cout << "# input-opt " << key.c_str() << " = " << val << std::endl;
}
////////////////////////////////////////////////////////////////////////
