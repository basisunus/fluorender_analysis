#include "TrackMap.h"
#include "tif_reader.h"
#include "lbl_reader.h"
#include <nrrd.h>
#include <iomanip>
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

int main(int argc, const char* argv[])
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
		if (tm_processor.ResolveForward(track_map, fi))
			printf("Frame %d and %d resolved.\n", fi, fi + 1);
		if (tm_processor.ResolveBackward(track_map, fi))
			printf("Frame %d and %d resolved.\n", fi, fi - 1);
	}

	printf("All done. Quitting.\n");

	return 0;
}