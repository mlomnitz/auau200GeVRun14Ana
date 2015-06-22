/* *********************************************************************
 *
 * Table of charm quark fragmentation ratios. 
 * Reference: ZEUS Collaboration - arXiv:hep-ex/0508019 - Table 4.
 *
 *  Authors:
 *            **Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  ** Code Maintainer
 *
 * *********************************************************************
*/
#include <iostream>
#include <utility>

using namespace std;

void printFragmentationRatios()
{
  cout<<"ZEUS Collaboration - arXiv:hep-ex/0508019 - Table 4."<<endl;
  cout<<"f(c → D+) 0.217 ± 0.014"<<endl;
  cout<<"f(c → D0) 0.523 ± 0.021"<<endl;
  cout<<"f(c → D+s) 0.095 ± 0.008"<<endl;
  cout<<"f(c → Λ+c) 0.144 ± 0.022"<<endl;
  cout<<"f(c → D∗+) 0.200 ± 0.009"<<endl;
}

namespace fragmentationRatios
{
  std::pair<double,double> const charmToDPlus          = make_pair(0.217,0.014); // value +- error
  std::pair<double,double> const charmToDZero          = make_pair(0.523,0.021);
  std::pair<double,double> const charmToDsPlus         = make_pair(0.095, 0.008);
  std::pair<double,double> const charmToLambdaCPlus    = make_pair(0.144,0.022);
  std::pair<double,double> const charmToDStarPlus      = make_pair(0.200,0.009);
}
