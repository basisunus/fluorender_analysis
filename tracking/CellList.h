#ifndef FL_CellList_h
#define FL_CellList_h

#include "Cell.h"
#include <boost\unordered_map.hpp>

namespace FL
{
	typedef boost::unordered_map<unsigned int, pCell> CellList;
	typedef boost::unordered_map<unsigned int, pCell>::iterator CellListIter;
}//namespace FL

#endif//FL_CellList_h