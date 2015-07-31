#include "TrackMap.h"

using namespace FL;

TrackMap::TrackMap() :
m_frame_num(0),
m_size_x(0),
m_size_y(0),
m_size_z(0),
m_data_bits(8),
m_scale(1.0f)
{
}

TrackMap::~TrackMap()
{
}

void TrackMapProcessor::SetSizes(TrackMap& track_map,
	size_t nx, size_t ny, size_t nz)
{
	track_map.m_size_x = nx;
	track_map.m_size_y = ny;
	track_map.m_size_z = nz;
}

void TrackMapProcessor::SetBits(TrackMap& track_map,
	size_t bits)
{
	track_map.m_data_bits = bits;
}

bool TrackMapProcessor::InitializeFrame(TrackMap& track_map,
	void* data, void* label, size_t frame)
{
	if (!data || !label)
		return false;

	//add one empty cell list to track_map
	track_map.m_cells_list.push_back(CellList());
	CellList &cell_list = track_map.m_cells_list.back();
	CellListIter iter;
	//in the meanwhile build the intra graph
	track_map.m_intra_graph_list.push_back(IntraGraph());

	size_t index;
	size_t i, j, k;
	size_t nx = track_map.m_size_x;
	size_t ny = track_map.m_size_y;
	size_t nz = track_map.m_size_z;
	float data_value;
	unsigned int label_value;

	//build cell list
	for (i = 0; i < nx; ++i)
	for (j = 0; j < ny; ++j)
	for (k = 0; k < nz; ++k)
	{
		index = nx*ny*k + nx*j + i;
		label_value = ((unsigned int*)label)[index];

		if (!label_value)
			continue;

		if (track_map.m_data_bits == 8)
			data_value = ((unsigned char*)data)[index] / 255.0f;
		else if (track_map.m_data_bits == 16)
			data_value = ((unsigned short*)data)[index] * track_map.m_scale / 65535.0f;

		iter = cell_list.find(label_value);
		if (iter != cell_list.end())
		{
			iter->second->Inc(i, j, k, data_value);
		}
		else
		{
			Cell* cell = new Cell(label_value);
			cell->Inc(i, j, k, data_value);
			cell_list.insert(std::pair<unsigned int, pCell>
				(label_value, pCell(cell)));
		}
	}

	//build intra graph
	for (i = 0; i < nx; ++i)
	for (j = 0; j < ny; ++j)
	for (k = 0; k < nz; ++k)
	{
		index = nx*ny*k + nx*j + i;
		label_value = ((unsigned int*)label)[index];
		if (!label_value)
			continue;

		iter = cell_list.find(label_value);
		if (iter != cell_list.end())
		{
			CheckCellContact(track_map,
				iter->second, data, label,
				i, j, k);
		}
	}

	//build vertex list
	track_map.m_vertices_list.push_back(VertexList());
	VertexList &vertex_list = track_map.m_vertices_list.back();
	for (CellList::iterator cell_iterator = cell_list.begin();
		cell_iterator != cell_list.end(); ++cell_iterator)
	{
		pVertex vertex(new Vertex(cell_iterator->second->Id()));
		vertex->SetCenter(cell_iterator->second->GetCenter());
		vertex->SetSize(cell_iterator->second->GetSizeUi(),
			cell_iterator->second->GetSizeF());
		vertex->AddCell(cell_iterator->second);
		cell_iterator->second->AddVertex(vertex);
		vertex_list.insert(std::pair<unsigned int, pVertex>
			(vertex->Id(), vertex));
	}

	track_map.m_frame_num++;
	return true;
}

#define ADD_CONTACT \
	idn = ((unsigned int*)label)[indexn]; \
	if (!idn) \
		ec++; \
		else if (idn != id) \
	{ \
		ec++; \
		if (data_bits == 8) \
			valuen = ((unsigned char*)data)[indexn] / 255.0f; \
				else if (data_bits == 16) \
			valuen = ((unsigned short*)data)[indexn] * scale / 65535.0f; \
		contact_value = std::min(value, valuen); \
		if (contact_value > m_contact_thresh) \
				{ \
			cc++; \
			iter = cell_list.find(idn); \
			if (iter != cell_list.end()) \
				AddContact(intra_graph, cell, iter->second, contact_value); \
		} \
	}

bool TrackMapProcessor::CheckCellContact(TrackMap& track_map,
	pCell &cell, void *data, void *label,
	size_t ci, size_t cj, size_t ck)
{
	int ec = 0;//external count
	int cc = 0;//contact count
	size_t nx = track_map.m_size_x;
	size_t ny = track_map.m_size_y;
	size_t nz = track_map.m_size_z;
	size_t indexn;//neighbor index
	unsigned int idn;//neighbor id
	float valuen;//neighbor vlaue
	size_t index = nx*ny*ck + nx*cj + ci;
	unsigned int id = cell->Id();
	float value;
	size_t data_bits = track_map.m_data_bits;
	float scale = track_map.m_scale;
	if (data_bits == 8)
		value = ((unsigned char*)data)[index] / 255.0f;
	else if (data_bits == 16)
		value = ((unsigned short*)data)[index] * scale / 65535.0f;
	float contact_value;
	IntraGraph &intra_graph = track_map.m_intra_graph_list.back();
	CellList &cell_list = track_map.m_cells_list.back();
	CellListIter iter;

	if (ci == 0)
		ec++;
	else
	{
		indexn = index - 1;
		ADD_CONTACT;
	}
	if (ci >= nx - 1)
		ec++;
	else
	{
		indexn = index + 1;
		ADD_CONTACT;
	}
	if (cj == 0)
		ec++;
	else
	{
		indexn = index - nx;
		ADD_CONTACT;
	}
	if (cj >= ny - 1)
		ec++;
	else
	{
		indexn = index + nx;
		ADD_CONTACT;
	}
	if (ck == 0)
		ec++;
	else
	{
		indexn = index - nx*ny;
		ADD_CONTACT;
	}
	if (ck >= nz - 1)
		ec++;
	else
	{
		indexn = index + nx*ny;
		ADD_CONTACT;
	}

	if (ec)
		cell->IncExternal(value);

	return true;
}

bool TrackMapProcessor::AddContact(IntraGraph& graph,
	pCell &cell1, pCell &cell2, float contact_value)
{
	IntraVert v1 = cell1->GetIntraVert();
	IntraVert v2 = cell2->GetIntraVert();
	if (v1 == IntraGraph::null_vertex())
	{
		v1 = boost::add_vertex(graph);
		graph[v1].id = cell1->Id();
		graph[v1].cell = cell1;
		cell1->SetIntraVert(v1);
	}
	if (v2 == IntraGraph::null_vertex())
	{
		v2 = boost::add_vertex(graph);
		graph[v2].id = cell2->Id();
		graph[v2].cell = cell2;
		cell2->SetIntraVert(v2);
	}

	std::pair<IntraEdge, bool> e = boost::edge(v1, v2, graph);
	if (!e.second)
	{
		e = boost::add_edge(v1, v2, graph);
		graph[e.first].size_ui = 1;
		graph[e.first].size_f = contact_value;
	}
	else
	{
		graph[e.first].size_ui++;
		graph[e.first].size_f += contact_value;
	}

	return true;
}

bool TrackMapProcessor::LinkMaps(TrackMap& track_map,
	size_t f1, size_t f2, void *data1, void *data2,
	void *label1, void *label2)
{
	size_t frame_num = track_map.m_frame_num;
	if (f1 >= frame_num || f2 >= frame_num ||
		f1 == f2 || !data1 || !data2 ||
		!label1 || !label2)
		return false;

	track_map.m_inter_graph_list.push_back(InterGraph());
	InterGraph &inter_graph = track_map.m_inter_graph_list.back();

	size_t index;
	size_t i, j, k;
	size_t nx = track_map.m_size_x;
	size_t ny = track_map.m_size_y;
	size_t nz = track_map.m_size_z;
	float data_value1, data_value2;
	unsigned int label_value1, label_value2;
	VertexList &vertex_list1 = track_map.m_vertices_list.at(f1);
	VertexList &vertex_list2 = track_map.m_vertices_list.at(f2);
	VertexListIter iter1, iter2;

	for (i = 0; i < nx; ++i)
	for (j = 0; j < ny; ++j)
	for (k = 0; k < nz; ++k)
	{
		index = nx*ny*k + nx*j + i;
		label_value1 = ((unsigned int*)label1)[index];
		label_value2 = ((unsigned int*)label2)[index];

		if (!label_value1 || !label_value2)
			continue;

		if (track_map.m_data_bits == 8)
		{
			data_value1 = ((unsigned char*)data1)[index] / 255.0f;
			data_value2 = ((unsigned char*)data2)[index] / 255.0f;
		}
		else if (track_map.m_data_bits == 16)
		{
			data_value1 = ((unsigned short*)data1)[index] * track_map.m_scale / 65535.0f;
			data_value2 = ((unsigned short*)data2)[index] * track_map.m_scale / 65535.0f;
		}

		iter1 = vertex_list1.find(label_value1);
		iter2 = vertex_list2.find(label_value2);

		if (iter1 == vertex_list1.end() ||
			iter2 == vertex_list2.end())
			continue;

		if (iter1->second->GetSizeF() < m_size_thresh ||
			iter2->second->GetSizeF() < m_size_thresh)
			continue;

		LinkVertices(inter_graph,
			iter1->second, iter2->second,
			std::min(data_value1, data_value2));
	}

	return true;
}

bool TrackMapProcessor::LinkVertices(InterGraph& graph,
	pVertex &vertex1, pVertex &vertex2, float overlap_value)
{
	InterVert v1 = vertex1->GetInterVert();
	InterVert v2 = vertex2->GetInterVert();
	if (v1 == InterGraph::null_vertex())
	{
		v1 = boost::add_vertex(graph);
		graph[v1].vertex = vertex1;
		vertex1->SetInterVert(v1);
	}
	if (v2 == InterGraph::null_vertex())
	{
		v2 = boost::add_vertex(graph);
		graph[v2].vertex = vertex2;
		vertex2->SetInterVert(v2);
	}

	std::pair<InterEdge, bool> e = boost::edge(v1, v2, graph);
	if (!e.second)
	{
		e = boost::add_edge(v1, v2, graph);
		graph[e.first].size_ui = 1;
		graph[e.first].size_f = overlap_value;
		FLIVR::Point p1 = vertex1->GetCenter();
		FLIVR::Point p2 = vertex2->GetCenter();
		graph[e.first].dist = float((p1 - p2).length());
		graph[e.first].link = 0;
	}
	else
	{
		graph[e.first].size_ui++;
		graph[e.first].size_f += overlap_value;
	}

	return true;
}

bool TrackMapProcessor::ResolveForward(TrackMap& track_map, size_t frame)
{
	if (frame >= track_map.m_frame_num - 1)
		return false;

	VertexList &vertex_list1 = track_map.m_vertices_list.at(frame);
	VertexList &vertex_list2 = track_map.m_vertices_list.at(frame + 1);
	CellList &cell_list2 = track_map.m_cells_list.at(frame + 1);
	InterGraph &inter_graph = track_map.m_inter_graph_list.at(frame);
	IntraGraph &intra_graph = track_map.m_intra_graph_list.at(frame + 1);

	VertexListIter iter;
	InterVert v1, v2;
	std::pair<InterIter, InterIter> adj_verts;
	InterIter inter_iter;
	std::vector<pwCell> cells;
	CellBinIter pwcell_iter;
	pVertex vertex2;
	std::vector<CellBin> cell_bins;
	pCell cell2, cell2c;
	IntraVert c2, c2c;
	std::pair<IntraIter, IntraIter> adj_cells;
	IntraIter intra_iter;
	bool added;
	std::pair<IntraEdge, bool> intra_edge;
	float osizef, c1sizef, c2sizef;
	size_t i;
	std::pair<InterEdge, bool> inter_edge, max_edge;
	float max_size_f;

	//check all vertices in the time frame
	for (iter = vertex_list1.begin();
		iter != vertex_list1.end(); ++iter)
	{
		cells.clear();
		v1 = iter->second->GetInterVert();
		adj_verts = boost::adjacent_vertices(v1, inter_graph);
		//for each adjacent vertex
		for (inter_iter = adj_verts.first;
			inter_iter != adj_verts.second; ++inter_iter)
		{
			v2 = *inter_iter;
			vertex2 = inter_graph[v2].vertex.lock();
			if (!vertex2)
				continue;
			//store all cells in the list temporarily
			for (pwcell_iter = vertex2->GetCellsBegin();
				pwcell_iter != vertex2->GetCellsEnd(); ++pwcell_iter)
				cells.push_back(*pwcell_iter);
		}
		//if a cell in the list has contacts that are also in the list,
		//try to group them
		for (pwcell_iter = cells.begin();
			pwcell_iter != cells.end(); ++pwcell_iter)
		{
			cell2 = pwcell_iter->lock();
			if (!cell2)
				continue;
			added = false;
			c2 = cell2->GetIntraVert();
			adj_cells = boost::adjacent_vertices(c2, intra_graph);
			//for each cell in contact
			for (intra_iter = adj_cells.first;
				intra_iter != adj_cells.second; ++intra_iter)
			{
				c2c = *intra_iter;
				if (FindCellBin(cells, intra_graph[c2c].cell))
				{
					intra_edge = boost::edge(c2, c2c, intra_graph);
					if (!intra_edge.second)
						continue;
					osizef = intra_graph[intra_edge.first].size_f;
					c1sizef = cell2->GetSizeF();
					cell2c = intra_graph[c2c].cell.lock();
					if (!cell2c)
						continue;
					c2sizef = cell2c->GetSizeF();
					if (osizef/c1sizef>m_contact_thresh ||
						osizef/c2sizef>m_contact_thresh)
					{
						//add both to bin list
						AddCellBin(cell_bins, *pwcell_iter, intra_graph[c2c].cell);
						added = true;
					}
				}
			}
			if (!added)//add to bin as well
				AddCellBin(cell_bins, *pwcell_iter);
		}

		//modify vertex list 2 if necessary
		for (i = 0; i < cell_bins.size(); ++i)
			MergeCells(vertex_list2, inter_graph, cell_bins[i]);

		//set flag for link
		adj_verts = boost::adjacent_vertices(v1, inter_graph);
		//for each adjacent vertex
		max_size_f = 0.0f;
		max_edge.second = false;
		for (inter_iter = adj_verts.first;
			inter_iter != adj_verts.second; ++inter_iter)
		{
			v2 = *inter_iter;
			inter_edge = boost::edge(v1, v2, inter_graph);
			if (inter_edge.second)
			{
				if (inter_graph[inter_edge.first].size_f > max_size_f)
				{
					max_size_f = inter_graph[inter_edge.first].size_f;
					max_edge = inter_edge;
				}
			}
		}
		if (max_edge.second)
			inter_graph[max_edge.first].link = 1;
	}

	return true;
}

bool TrackMapProcessor::ResolveBackward(TrackMap& track_map, size_t frame)
{
	return true;
}

bool TrackMapProcessor::EqualCells(pwCell &cell1, pwCell &cell2)
{
	return !cell1.owner_before(cell2) && !cell2.owner_before(cell1);
}

bool TrackMapProcessor::FindCellBin(CellBin &bin, pwCell &cell)
{
	for (size_t i = 0; i < bin.size(); ++i)
		if (EqualCells(bin[i], cell))
			return true;
	return false;
}

bool TrackMapProcessor::AddCellBin(std::vector<CellBin> &bins, pwCell &cell)
{
	bool found_cell;
	for (size_t i = 0; i < bins.size(); ++i)
	{
		found_cell = FindCellBin(bins[i], cell);
		if (found_cell)
			return false;
	}
	CellBin bin;
	bin.push_back(cell);
	bins.push_back(bin);
	return true;
}

bool TrackMapProcessor::AddCellBin(std::vector<CellBin> &bins, pwCell &cell1, pwCell &cell2)
{
	bool found_cell1, found_cell2;
	for (size_t i = 0; i < bins.size(); ++i)
	{
		found_cell1 = FindCellBin(bins[i], cell1);
		found_cell2 = FindCellBin(bins[i], cell2);
		if (found_cell1 && found_cell2)
			return false;
		else if (found_cell1 && !found_cell2)
		{
			bins[i].push_back(cell2);
			return true;
		}
		else if (!found_cell1 && found_cell2)
		{
			bins[i].push_back(cell1);
			return true;
		}
	}
	CellBin bin;
	bin.push_back(cell1);
	bin.push_back(cell2);
	bins.push_back(bin);
	return true;
}

size_t TrackMapProcessor::GetBinsCellCount(std::vector<CellBin> &bins)
{
	size_t count = 0;
	for (size_t i = 0; i < bins.size(); ++i)
		count += bins[i].size();
	return count;
}

bool TrackMapProcessor::MergeCells(VertexList& vertex_list,
	InterGraph &graph, CellBin &bin)
{
	if (bin.size() <= 1)
		return false;

	//the keeper
	pCell cell0, cell;
	pVertex vertex0, vertex;
	VertexListIter vert_iter;
	std::pair<InterIter, InterIter> adj_verts;
	InterVert inter_vert;
	InterIter inter_iter;
	std::pair<InterEdge, bool> e, e0;

	for (size_t i = 0; i < bin.size(); ++i)
	{
		cell = bin[i].lock();
		if (!cell)
			continue;
		if (!vertex0)
		{
			cell0 = cell;
			vertex0 = cell0->GetVertex().lock();
		}
		else
		{
			vertex = cell->GetVertex().lock();
			if (!vertex)
				continue;
			//relink inter graph
			inter_vert = vertex->GetInterVert();
			adj_verts = boost::adjacent_vertices(inter_vert, graph);
			//for each adjacent vertex
			for (inter_iter = adj_verts.first;
				inter_iter != adj_verts.second; ++inter_iter)
			{
				//get edge
				e = boost::edge(inter_vert, *inter_iter, graph);
				if (!e.second)
					continue;
				//add an edge between vertex0 and inter_iter
				e0 = boost::edge(*inter_iter,
					vertex0->GetInterVert(), graph);
				if (!e0.second)
				{
					e0 = boost::add_edge(*inter_iter,
						vertex0->GetInterVert(), graph);
					graph[e0.first].size_ui = graph[e.first].size_ui;
					graph[e0.first].size_f = graph[e.first].size_f;
					graph[e0.first].dist = graph[e.first].dist;
					graph[e0.first].link = graph[e.first].link;
				}
				else
				{
					graph[e0.first].size_ui += graph[e.first].size_ui;
					graph[e0.first].size_f += graph[e.first].size_f;
				}
				//delete the old edge
				//graph.remove_edge(e.first);
			}
			//remove the vertex from inter graph
			//edges should be removed as well
			boost::remove_vertex(inter_vert, graph);

			//remove vertex from list
			vert_iter = vertex_list.find(vertex->Id());
			if (vert_iter != vertex_list.end())
				vertex_list.erase(vert_iter);

			//add cell to vertex0
			vertex0->AddCell(cell);
			cell->AddVertex(vertex0);
		}
	}

	return true;
}