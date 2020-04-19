////////////////////////////////////////////////////////////////////////
//#include <iostream>
//#include <sstream>
//#include <cmath>
#include <cstdio>
//#include <cstdlib>
//#include <stdio.h>
#include <ctime>
//#include <iterator>
//#include <algorithm>
////////////////////////////////////////////////////////////////////////
void print_date(){
  time_t now = time(NULL);
  struct tm *pnow = localtime(&now);

  char buff[128]="";

  printf("%2d:%2d:%2d:%2d:%2d:%2d\n", 
	 pnow->tm_year+199, 
	 pnow->tm_mon+1, 
	 pnow->tm_mday, 
	 pnow->tm_hour, 
	 pnow->tm_min, 
	 pnow->tm_sec);

//  sprintf(buff,"%d:%d:%d",pnow->tm_hour,pnow->tm_min,pnow->tm_sec);
//  printf("\n########\n");
//  printf(buff);
//  printf("\n########\n");
//  getchar();
//  namespace pt = boost::posix_time;
//  namespace gg = boost::gregorian;
//
//  typedef boost::date_time::c_local_adjustor<pt::ptime> local_adj;
//  auto epoch = local_adj::utc_to_local(pt::ptime(gg::date(1970, 1, 1)));
//
//  auto facet = new pt::time_facet("%Y-%m-%d %H:%M:%S");
//  std::stringstream ss;
//  ss.imbue(std::locale(std::cout.getloc(), facet));
//        
//  pt::seconds time_stamp(1371292653L);
//  auto date_time = epoch + time_stamp;
//  ss << date_time;
//  std::cout << ss.str() << std::endl;
//
//  return 0;
}
////////////////////////////////////////////////////////////////////////
