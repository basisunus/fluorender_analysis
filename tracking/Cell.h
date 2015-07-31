#ifndef FL_Cell_h
#define FL_Cell_h

#include <Point.h>
#include <boost\shared_ptr.hpp>
#include <boost\weak_ptr.hpp>
#include <boost\graph\graph_traits.hpp>
#include <boost\graph\adjacency_list.hpp>
#include <boost\graph\adjacency_iterator.hpp>

namespace FL
{
	class Vertex;
	typedef boost::shared_ptr<Vertex> pVertex;
	typedef boost::weak_ptr<Vertex> pwVertex;
	class Cell;
	typedef boost::shared_ptr<Cell> pCell;
	typedef boost::weak_ptr<Cell> pwCell;

	struct IntraEdgeData
	{
		unsigned int size_ui;
		float size_f;
	};

	struct IntraCellData
	{
		unsigned int id;
		pwCell cell;
	};

	typedef boost::adjacency_list<boost::vecS,
		boost::vecS, boost::undirectedS,
		IntraCellData, IntraEdgeData> IntraGraph;
	typedef IntraGraph::vertex_descriptor IntraVert;
	typedef IntraGraph::edge_descriptor IntraEdge;
	typedef boost::graph_traits<IntraGraph>::adjacency_iterator IntraIter;

	class Cell
	{
	public:
		Cell(unsigned int id) :
			m_id(id), m_size_ui(0), m_size_f(0.0f),
			m_external_ui(0), m_external_f(0.0f),
			m_intra_vert(IntraGraph::null_vertex())
		{}
		~Cell() {}

		unsigned int Id();
		IntraVert GetIntraVert();
		void SetIntraVert(IntraVert intra_vert);
		void Inc(size_t i, size_t j, size_t k, float value);
		void IncExternal(float value);
		void AddVertex(pVertex &vertex);
		pwVertex GetVertex();

		FLIVR::Point &GetCenter();
		unsigned int GetSizeUi();
		float GetSizeF();

	private:
		unsigned int m_id;
		FLIVR::Point m_center;
		//size
		unsigned int m_size_ui;
		float m_size_f;
		//external size
		unsigned int m_external_ui;
		float m_external_f;
		IntraVert m_intra_vert;
		pwVertex m_vertex;//parent
	};

	inline unsigned int Cell::Id()
	{
		return m_id;
	}

	inline IntraVert Cell::GetIntraVert()
	{
		return m_intra_vert;
	}

	inline void Cell::SetIntraVert(IntraVert intra_vert)
	{
		m_intra_vert = intra_vert;
	}

	inline void Cell::Inc(size_t i, size_t j, size_t k, float value)
	{
		m_center = FLIVR::Point(
			(m_center*m_size_ui + FLIVR::Point(double(i),
			double(j), double(k))) / (m_size_ui + 1));
		m_size_ui++;
		m_size_f += value;
	}

	inline void Cell::IncExternal(float value)
	{
		m_external_ui++;
		m_external_f += value;
	}

	inline void Cell::AddVertex(pVertex &vertex)
	{
		m_vertex = vertex;
	}

	inline pwVertex Cell::GetVertex()
	{
		return m_vertex;
	}

	inline FLIVR::Point &Cell::GetCenter()
	{
		return m_center;
	}

	inline unsigned int Cell::GetSizeUi()
	{
		return m_size_ui;
	}

	inline float Cell::GetSizeF()
	{
		return m_size_f;
	}

}//namespace FL

#endif//FL_Cell_h