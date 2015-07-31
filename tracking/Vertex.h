#ifndef FL_Vertex_h
#define FL_Vertex_h

#include "Cell.h"
#include <vector>

namespace FL
{
	struct InterEdgeData
	{
		unsigned int size_ui;
		float size_f;
		float dist;
		unsigned int link;
	};

	struct InterVertexData
	{
		pwVertex vertex;
	};

	typedef boost::adjacency_list<boost::vecS,
		boost::vecS, boost::undirectedS,
		InterVertexData, InterEdgeData> InterGraph;
	typedef InterGraph::vertex_descriptor InterVert;
	typedef InterGraph::edge_descriptor InterEdge;
	typedef boost::graph_traits<InterGraph>::adjacency_iterator InterIter;

	typedef std::vector<pwCell> CellBin;
	typedef CellBin::iterator CellBinIter;

	class Vertex
	{
	public:
		Vertex(unsigned int id) :
			m_id(id), m_size_ui(0), m_size_f(0.0f),
			m_inter_vert(InterGraph::null_vertex())
		{}
		~Vertex() {};

		unsigned int Id();
		InterVert GetInterVert();
		void SetInterVert(InterVert inter_vert);

		void SetCenter(FLIVR::Point &center);
		void SetSize(unsigned int size_ui, float size_f);

		//cells
		void AddCell(pCell &cell);
		CellBinIter GetCellsBegin();
		CellBinIter GetCellsEnd();

		FLIVR::Point &GetCenter();
		unsigned int GetSizeUi();
		float GetSizeF();

	private:
		unsigned int m_id;
		FLIVR::Point m_center;
		unsigned int m_size_ui;
		float m_size_f;
		InterVert m_inter_vert;
		CellBin m_cells;//children
	};

	inline unsigned int Vertex::Id()
	{
		return m_id;
	}

	inline InterVert Vertex::GetInterVert()
	{
		return m_inter_vert;
	}

	inline void Vertex::SetInterVert(InterVert inter_vert)
	{
		m_inter_vert = inter_vert;
	}

	inline void Vertex::SetCenter(FLIVR::Point &center)
	{
		m_center = center;
	}

	inline void Vertex::SetSize(unsigned int size_ui, float size_f)
	{
		m_size_ui = size_ui;
		m_size_f = size_f;
	}

	inline void Vertex::AddCell(pCell &cell)
	{
		m_cells.push_back(pwCell(cell));
	}

	inline CellBinIter Vertex::GetCellsBegin()
	{
		return m_cells.begin();
	}

	inline CellBinIter Vertex::GetCellsEnd()
	{
		return m_cells.end();
	}

	inline FLIVR::Point &Vertex::GetCenter()
	{
		return m_center;
	}

	inline unsigned int Vertex::GetSizeUi()
	{
		return m_size_ui;
	}

	inline float Vertex::GetSizeF()
	{
		return m_size_f;
	}

}//namespace FL

#endif//FL_Vertex_h