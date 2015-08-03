#ifndef FL_TrackMap_h
#define FL_TrackMap_h

#include "CellList.h"
#include "VertexList.h"

namespace FL
{
//tags
#define TAG_CELL		1
#define TAG_EDGE		2
#define TAG_VERT		3
#define TAG_CMAP		4
#define TAG_FRAM		5

	class TrackMap;
	class TrackMapProcessor
	{
	public:
		TrackMapProcessor() :
		m_contact_thresh(0.2f), m_size_thresh(30.0f) {};
		~TrackMapProcessor() {};

		void SetContactThresh(float value);
		void SetSizeThresh(float value);

		void SetSizes(TrackMap& track_map,
			size_t nx, size_t ny, size_t nz);
		void SetBits(TrackMap& track_map,
			size_t bits);

		bool InitializeFrame(TrackMap& track_map,
			void *data, void *label, size_t frame);
		bool LinkMaps(TrackMap& track_map,
			size_t f1, size_t f2,
			void *data1, void *data2,
			void *label1, void *label2);
		bool ResolveForward(TrackMap& track_map, size_t frame);
		bool ResolveBackward(TrackMap& track_map, size_t frame);

		bool Export(TrackMap& track_map, std::string &filename);

	private:
		float m_contact_thresh;
		float m_size_thresh;

		bool CheckCellContact(TrackMap& track_map,
			pCell &cell, void *data, void *label,
			size_t ci, size_t cj, size_t ck);
		bool AddContact(IntraGraph& graph,
			pCell &cell1, pCell &cell2,
			float contact_value);
		bool LinkVertices(InterGraph& graph,
			pVertex &vertex1, pVertex &vertex2,
			float overlap_value);
		bool EqualCells(pwCell &cell1, pwCell &cell2);
		bool FindCellBin(CellBin &bin, pwCell &cell);
		bool AddCellBin(std::vector<CellBin> &bins,
			pwCell &cell);
		bool AddCellBin(std::vector<CellBin> &bins,
			pwCell &cell1, pwCell &cell2);
		size_t GetBinsCellCount(std::vector<CellBin> &bins);
		bool MergeCells(VertexList& vertex_list,
			InterGraph &graph, CellBin &bin);
	};

	inline void TrackMapProcessor::SetContactThresh(float value)
	{
		m_contact_thresh = value;
	}

	inline void TrackMapProcessor::SetSizeThresh(float value)
	{
		m_size_thresh = value;
	}

	class TrackMap
	{
	public:
		TrackMap();
		~TrackMap();

		size_t GetFrameNum();

	private:
		//data information
		size_t m_frame_num;
		size_t m_size_x;
		size_t m_size_y;
		size_t m_size_z;
		size_t m_data_bits;
		float m_scale;

		//lists
		std::vector<CellList> m_cells_list;
		std::vector<VertexList> m_vertices_list;
		std::vector<IntraGraph> m_intra_graph_list;
		std::vector<InterGraph> m_inter_graph_list;

		friend class TrackMapProcessor;
	};

	inline size_t TrackMap::GetFrameNum()
	{
		return m_frame_num;
	}

}//namespace FL

#endif//FL_TrackMap_h