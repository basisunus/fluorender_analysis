#include "TrackMap.h"
#include "tif_reader.h"
#include "lbl_reader.h"
#include "msk_writer.h"
#include <nrrd.h>
#include <iomanip>
#include <codecvt>
#include <boost/unordered_map.hpp>
#include <vld.h>

using namespace std;

Nrrd* ReadData(string &filename)
{
	TIFReader reader;

	reader.SetFile(filename);
	reader.SetSliceSeq(false);
	wstring time_id = L"NA!@";
	reader.SetTimeId(time_id);
	reader.Preprocess();
	return reader.Convert(0, 0, true);
}

Nrrd* ReadLabel(string &filename)
{
	LBLReader reader;

	reader.SetFile(filename);
	return reader.Convert(0, 0, true);
}

inline wstring s2ws(const string& utf8)
{
	std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>, wchar_t> converter;
	return converter.from_bytes(utf8);
}

void WriteLabel(Nrrd* label, string &filename)
{
	MSKWriter writer;
	writer.SetData(label);
	writer.Save(s2ws(filename), 1);
}

unsigned int GetMappedID(unsigned int id, unsigned int* data_label1,
	unsigned int* data_label2, size_t size)
{
	for (size_t i = 0; i < size; ++i)
		if (data_label1[i] == id)
			return data_label2[i];
	return 0;
}

typedef boost::unordered_map<unsigned int, unsigned int> CellMap;
typedef boost::unordered_map<unsigned int, unsigned int>::iterator CellMapIter;

int main(int argc, const char* argv[])
{
	if (argc != 7)
	{
		printf("Wrong arguments.\n");
		return 0;
	}

	string in_format = argv[1];
	string out_format = argv[2];
	string mapfile = argv[3];
	int fdigits = atoi(argv[4]);
	int fstart = atoi(argv[5]);
	int fend = atoi(argv[6]);

	FL::TrackMap track_map;
	FL::TrackMapProcessor tm_processor;

	if (!tm_processor.Import(track_map, mapfile))
	{
		printf("Track map import error.\n");
		return 0;
	}

	ostringstream oss;
	string fn_data;
	string fn_label;
	string str;

	Nrrd* nrrd_label_in1 = 0;
	Nrrd* nrrd_label_in2 = 0;
	Nrrd* nrrd_label_out1 = 0;
	Nrrd* nrrd_label_out2 = 0;

	unsigned int *data_in1;
	unsigned int *data_in2;
	unsigned int *data_out1;
	unsigned int *data_out2;

	//read first frame
	oss.str("");
	oss << setfill('0') << setw(fdigits) << fstart;
	str = oss.str();
	nrrd_label_in1 = ReadLabel(in_format + str + ".lbl");
	if (!nrrd_label_in1)
	{
		printf("File error.\n");
		return 0;
	}
	size_t nx = nrrd_label_in1->axis[0].size;
	size_t ny = nrrd_label_in1->axis[1].size;
	size_t nz = nrrd_label_in1->axis[2].size;
	size_t size = nx*ny*nz;
	size_t i;

	CellMap cell_map;
	CellMapIter iter;
	unsigned int label_in1, label_in2;
	//duplicate
	nrrdCopy(nrrd_label_in2, nrrd_label_in1);
	for (i = 0; i < size; ++i)
	{
		label_in1 = ((unsigned int*)(nrrd_label_in1->data))[i];
		iter = cell_map.find(label_in1);
		if (iter != cell_map.end())
		{
			((unsigned int*)(nrrd_label_in2->data))[i] = iter->second;
		}
		else
		{
			if (tm_processor.GetMappedID(track_map,
				label_in1, label_in2, 0))
			{
				((unsigned int*)(nrrd_label_in2->data))[i] = label_in2;
				cell_map.insert(pair<unsigned int, unsigned int>(label_in1, label_in2));
			}
		}
	}
	WriteLabel(nrrd_label_in2, out_format + str + ".lbl");

	unsigned int label_out1, label_out2;
	//remaining frames
	for (int fi = fstart + 1; fi <= fend; ++fi)
	{
		cell_map.clear();

		oss.str("");
		oss << setfill('0') << setw(fdigits) << fi;
		str = oss.str();
		nrrd_label_out1 = ReadLabel(in_format + str + ".lbl");
		nrrdCopy(nrrd_label_out2, nrrd_label_out1);

		for (i = 0; i < size; ++i)
		{
			label_out1 = ((unsigned int*)(nrrd_label_out1->data))[i];
			iter = cell_map.find(label_out1);
			if (iter != cell_map.end())
			{
				((unsigned int*)(nrrd_label_out2->data))[i] = iter->second;
			}
			else
			{
				if (tm_processor.GetMappedID(track_map,
					label_out1, label_in1, fi, fi - 1))
				{
					label_out2 = GetMappedID(label_in1,
						(unsigned int*)(nrrd_label_in1->data),
						(unsigned int*)(nrrd_label_in2->data),
						size);
					if (label_out2)
					{
						((unsigned int*)(nrrd_label_out2->data))[i] = label_out2;
						cell_map.insert(pair<unsigned int, unsigned int>(label_out1, label_out2));
					}
				}
			}
		}

		//save
		WriteLabel(nrrd_label_out2, out_format + str + ".lbl");

		//swap
		nrrdNuke(nrrd_label_in1);
		nrrdNuke(nrrd_label_in2);
		nrrd_label_in1 = nrrd_label_out1;
		nrrd_label_in2 = nrrd_label_out2;
	}

	nrrdNuke(nrrd_label_out1);
	nrrdNuke(nrrd_label_out2);

	return 0;
}

/*int main(int argc, const char* argv[])
{
	if (argc != 6)
	{
		printf("Wrong arguments.\n");
		return 0;
	}

	string format = argv[1];
	string outfilename = argv[2];
	int fdigits = atoi(argv[3]);
	int fstart = atoi(argv[4]);
	int fend = atoi(argv[5]);

	FL::TrackMap track_map;
	FL::TrackMapProcessor tm_processor;

	ostringstream oss;
	string fn_data;
	string fn_label;
	string str;

	Nrrd* nrrd_data1 = 0;
	Nrrd* nrrd_data2 = 0;
	Nrrd* nrrd_label1 = 0;
	Nrrd* nrrd_label2 = 0;

	//read first frame
	oss.str("");
	oss << setfill('0') << setw(fdigits) << fstart;
	str = oss.str();
	nrrd_data1 = ReadData(format + str + ".tif");
	nrrd_label1 = ReadLabel(format + str + ".lbl");
	if (!nrrd_data1 || !nrrd_label1)
	{
		printf("File error.\n");
		return 0;
	}
	tm_processor.SetSizes(track_map,
		nrrd_data1->axis[0].size,
		nrrd_data1->axis[1].size,
		nrrd_data1->axis[2].size);
	//tm_processor.SetBits(track_map,
	//	8);
	tm_processor.SetContactThresh(0.2f);
	if (tm_processor.InitializeFrame(track_map,
		nrrd_data1->data, nrrd_label1->data, fstart))
		printf("Frame %d initialized.\n", fstart);
	else
		return 0;


	//initialization
	for (int fi = fstart+1; fi <= fend; ++fi)
	{
		oss.str("");
		oss << setfill('0') << setw(fdigits) << fi;
		str = oss.str();
		nrrd_data2 = ReadData(format+str+".tif");
		nrrd_label2 = ReadLabel(format+str+".lbl");

		if (!nrrd_data2 || !nrrd_label2)
		{
			printf("File error.\n");
			return 0;
		}

		if (tm_processor.InitializeFrame(track_map,
			nrrd_data2->data, nrrd_label2->data, fi))
			printf("Frame %d initialized.\n", fi);
		else
			return 0;

		//link maps 1 and 2
		if (tm_processor.LinkMaps(track_map, fi - 1, fi,
			nrrd_data1->data, nrrd_data2->data,
			nrrd_label1->data, nrrd_label2->data))
			printf("Frame %d and %d linked.\n", fi - 1, fi);
		else
			return 0;

		nrrdNuke(nrrd_data1);
		nrrdNuke(nrrd_label1);
		nrrd_data1 = nrrd_data2;
		nrrd_label1 = nrrd_label2;
	}

	nrrdNuke(nrrd_data2);
	nrrdNuke(nrrd_label2);

	//resolve multiple links of single vertex
	for (size_t fi = 0; fi < track_map.GetFrameNum(); ++fi)
	{
		if (tm_processor.ResolveGraph(track_map, fi, fi + 1))
			printf("Frame %zd and %zd resolved.\n", fi, fi + 1);
		if (tm_processor.ResolveGraph(track_map, fi, fi - 1))
			printf("Frame %zd and %zd resolved.\n", fi, fi - 1);
	}

	//save
	if (tm_processor.Export(track_map, outfilename))
		printf("Track map saved.\n");

	printf("All done. Quitting.\n");

	return 0;
}*/

