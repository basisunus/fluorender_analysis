#ifndef FL_VertexList_h
#define FL_VertexList_h

#include "Vertex.h"
#include <boost\unordered_map.hpp>

namespace FL
{
	typedef boost::unordered_map<unsigned int, pVertex> VertexList;
	typedef boost::unordered_map<unsigned int, pVertex>::iterator VertexListIter;
}//namespace FL

#endif//FL_VertexList_h