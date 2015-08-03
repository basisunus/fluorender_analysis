#ifndef FL_TrackMap_h
#define FL_TrackMap_h

#include "CellList.h"
#include "VertexList.h"
#include <fstream>

namespace FL
{
//tags
#define TAG_CELL		1
#define TAG_VERT		2
#define TAG_INTRA_EDGE	3
#define TAG_INTER_EDGE	4
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
		bool ResolveGraph(TrackMap& track_map, size_t frame1, size_t frame2);

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

		//export
		void WriteBool(std::ofstream& ofs, bool value);
		void WriteTag(std::ofstream& ofs, unsigned char tag);
		void WriteUint(std::ofstream& ofs, unsigned int value);
		void WriteFloat(std::ofstream& ofs, float value);
		void WritePoint(std::ofstream& ofs, FLIVR::Point &point);
		void WriteCell(std::ofstream& ofs, pCell &cell);
		void WriteVertex(std::ofstream& ofs, pVertex &vertex);
	};

	inline void TrackMapProcessor::SetContactThresh(float value)
	{
		m_contact_thresh = value;
	}

	inline void TrackMapProcessor::SetSizeThresh(float value)
	{
		m_size_thresh = value;
	}

	inline void TrackMapProcessor::WriteBool(std::ofstream& ofs, bool value)
	{
		ofs.write(reinterpret_cast<const char*>(&value), sizeof(bool));
	}

	inline void TrackMapProcessor::WriteTag(std::ofstream& ofs, unsigned char tag)
	{
		ofs.write(reinterpret_cast<const char*>(&tag), sizeof(unsigned char));
	}

	inline void TrackMapProcessor::WriteUint(std::ofstream& ofs, unsigned int value)
	{
		ofs.write(reinterpret_cast<const char*>(&value), sizeof(unsigned int));
	}

	inline void TrackMapProcessor::WriteFloat(std::ofstream& ofs, float value)
	{
		ofs.write(reinterpret_cast<const char*>(&value), sizeof(float));
	}

	inline void TrackMapProcessor::WritePoint(std::ofstream& ofs, FLIVR::Point &point)
	{
		double x = point.x();
		ofs.write(reinterpret_cast<const char*>(&x), sizeof(double));
		x = point.y();
		ofs.write(reinterpret_cast<const char*>(&x), sizeof(double));
		x = point.z();
		ofs.write(reinterpret_cast<const char*>(&x), sizeof(double));
	}

	inline void TrackMapProcessor::WriteCell(std::ofstream& ofs, pCell &cell)
	{
		WriteTag(ofs, TAG_CELL);
		WriteUint(ofs, cell->Id());
		WriteUint(ofs, cell->GetSizeUi());
		WriteFloat(ofs, cell->GetSizeF());
		WriteUint(ofs, cell->GetExternalUi());
		WriteFloat(ofs, cell->GetExternalF());
		WritePoint(ofs, cell->GetCenter());
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