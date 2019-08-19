#include "stdafx.h"

#include <iostream>
#include <fstream>

#include <vector>
#include <cmath>

#include <sstream>
#include <algorithm>
#include <ctime>

#include <omp.h>

#include "PSRFM.h"


using namespace Eigen;
using namespace std;


PSRFM::PSRFM()
{
	// do nothing
}


PSRFM::~PSRFM()
{
}

// delimiter: by it to split the input string
// position: 1, 2, 3....if known. 1000 the last splitted word (string)
string PSRFM::string_split(string input_str, char delimiter, int position) {

	string word;

	vector<string> words;

	istringstream iss(input_str);

	int words_count;
	int get_pos;

	while (getline(iss, word, delimiter)) {
		words.push_back(word);
	}

	words_count = static_cast<int>(words.size());

	get_pos = position - 1;
	if (position >= words_count) {
		get_pos = words_count - 1;
	}
	
	return words[get_pos];
}

string PSRFM::get_short_fine_sensor_name(string input_str) {
	
	char delimiter;
	string full_sensor_name;
	string short_sensor_name;


	// get rid of directories
	delimiter = '\\';
	full_sensor_name = string_split(input_str, delimiter, 1000);

	// get the first part of the file name
	delimiter = '_';
	short_sensor_name = string_split(full_sensor_name, delimiter, 1);

	return short_sensor_name;
}



void PSRFM::clean_memory(SENSOR_PAIR *sensor_p[], int PAIRS, CONTROL_PARAMETER *ctrl, SENSOR *pimage) {

	for (int p = 0; p < PAIRS; p++) {
		for (int b = 0; b < ctrl->NUM_BANDS; b++) {
			for (int r = 0; r < ctrl->NUM_ROWS; r++) {
				delete[] sensor_p[p]->fimage.data[b][r];
				delete[] sensor_p[p]->cimage.data[b][r];
			}
			delete[] sensor_p[p]->fimage.data[b];
			delete[] sensor_p[p]->cimage.data[b];
		}
		delete[] sensor_p[p]->fimage.data;
		delete[] sensor_p[p]->cimage.data;
	}
	for (int p = 0; p < PAIRS; p++) {
		for (int r = 0; r < ctrl->NUM_ROWS; r++) {
			delete[] sensor_p[p]->fimage.mask[r];
		}
		delete[] sensor_p[p]->fimage.mask;
	}

	delete ctrl;
	delete[] pimage;
	for (int p = 0; p < PAIRS; p++) {
		delete sensor_p[p];
	}
}


double PSRFM::correlation(float **x1, float **x2, int row, int col, int inval) {

	int i, j;
	int count;
	double mean_x1;
	double mean_x2;

	double x1mx2;

	double x1_sq;
	double x2_sq;
	double std_x1;
	double std_x2;

	double cc;

	count = 1;
	mean_x1 = 0.0;
	mean_x2 = 0.0;
	x1_sq = 0.0;
	x2_sq = 0.0;
	x1mx2 = 0.0;


	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			// deal with invalid pixels
			if (x1[i][j] != inval && x2[i][j] != inval) {
				mean_x1 += (x1[i][j] - mean_x1) / count;
				mean_x2 += (x2[i][j] - mean_x2) / count;
				x1mx2 += x1[i][j] * x2[i][j];
				x1_sq += x1[i][j] * x1[i][j];
				x2_sq += x2[i][j] * x2[i][j];
				count++;
			}
		}
	}

	std_x1 = sqrt(x1_sq - (count - 1) * mean_x1 * mean_x1);
	std_x2 = sqrt(x2_sq - (count - 1) * mean_x2 * mean_x2);

	cc = x1mx2 - (count - 1) * mean_x1 * mean_x2;
	cc /= (std_x1 * std_x2);

	return cc;
}

void PSRFM::initialize_ctrl_parameters(CONTROL_PARAMETER *ctrl) {
	ctrl->NUM_ROWS = 0;
	ctrl->NUM_COLS = 0;
	ctrl->RESOLUTION = 0;
	ctrl->NUM_PAIRS = 2;
	ctrl->CLUSTER_RANGE[0] = 8;
	ctrl->CLUSTER_RANGE[1] = 14;
	ctrl->CLUSTER_METHOD = "KMEAN";
	ctrl->CLUSTER_DATA = "fine";
	ctrl->CLUSTER_OPTIMAL = "CF";
	ctrl->MERGE_METHOD = "temporal";
}

void PSRFM::initialize_sensors(SENSOR_PAIR *sensor_p[], int PAIRS) {

	for (int i = 0; i < PAIRS; i++) {

		//fine
		sensor_p[i]->fimage.fname = "";
		sensor_p[i]->fimage.mfname = "";
		sensor_p[i]->fimage.fillv = -9999;
		sensor_p[i]->fimage.range[0] = 0;
		sensor_p[i]->fimage.range[1] = 10000;
		sensor_p[i]->fimage.res = 0;
		sensor_p[i]->fimage.scale = 10000;
		sensor_p[i]->fimage.uncertainty = 10000;

		//coarse 
		sensor_p[i]->cimage.fname = "";
		sensor_p[i]->cimage.mfname = "";
		sensor_p[i]->cimage.fillv = -9999;
		sensor_p[i]->cimage.range[0] = 0;
		sensor_p[i]->cimage.range[1] = 10000;
		sensor_p[i]->cimage.res = 0;
		sensor_p[i]->cimage.scale = 10000;
		sensor_p[i]->cimage.uncertainty = 10000;
	}
}

int PSRFM::read_image_data(SENSOR_PAIR *sensor_p[], CONTROL_PARAMETER *ctrl, int PAIRS, ofstream& log_file) {

	int i;
	int j;
	int b;
	int p;
	long inval_fpc = 0;
	long inval_cpc = 0;

	short ***cdata;

	string fname;
	ifstream fdataFStream;  // fine data file
	ifstream cdataFStream;	// coarse data file
	ifstream mfdataFStream;	// fine mask data file

	int num_bands;
	int bandSize = ctrl->NUM_ROWS*ctrl->NUM_COLS;
	int coarse_bandSize = ctrl->COARSE_ROWS*ctrl->COARSE_COLS;

	int image_date_int = 0;
	string image_date_str = "";

	// allocate memory for co-registration
	cdata = new short**[ctrl->NUM_BANDS];
	for (int b = 0; b < ctrl->NUM_BANDS; b++) {
		cdata[b] = new short*[ctrl->COARSE_ROWS];
		for (int r = 0; r < ctrl->COARSE_ROWS; r++) {
			cdata[b][r] = new short[ctrl->COARSE_COLS];
		}
	}

	for (p = 0; p < PAIRS; p++) {

		// fine image data
		fname = sensor_p[p]->fimage.fname;

		fdataFStream.open(fname.c_str(), ios::in | ios::binary | ios::ate);

		if (!fdataFStream.is_open()) {
			cout << "Cannot open input image file " << fname.c_str() << endl;
			return 1;
		}
		else { // read data
			num_bands = (int)fdataFStream.tellg() * sizeof(char) / (2 * bandSize);
			if (num_bands != ctrl->NUM_BANDS) {
				cout << "File not in correct dimension <" << fname.c_str() << ">" << endl;
			}
			else {
				// move to start of input File
				fdataFStream.seekg(0, ios_base::beg);
			}
		}

		// coarse image data
		fname = sensor_p[p]->cimage.fname;
		cdataFStream.open(fname.c_str(), ios::in | ios::binary | ios::ate);

		if (!cdataFStream.is_open()) {
			cout << "Cannot open input image file " << fname.c_str() << endl;
			return 1;
		}
		else { // read data
			num_bands = (int)cdataFStream.tellg() * sizeof(char) / (2 * coarse_bandSize);
			if (num_bands != ctrl->NUM_BANDS) {
				cout << "File not in correct dimension <" << fname.c_str() << ">" << endl;
			}
			else {
				// move to start of input File
				cdataFStream.seekg(0, ios_base::beg);
			}
		}

		// fine mask image data
		fname = sensor_p[p]->fimage.mfname;
		mfdataFStream.open(fname.c_str(), ios::in | ios::binary | ios::ate);

		if (!mfdataFStream.is_open()) {
			cout << "Cannot open input image file " << fname.c_str() << endl;
			return 1;
		}
		else { // read data
			num_bands = (int)mfdataFStream.tellg() * sizeof(char) / bandSize;
			if (num_bands != 1) {
				cout << "File not in correct dimension <" << fname.c_str() << ">" << endl;
			}
			else {
				// move to start of input File
				mfdataFStream.seekg(0, ios_base::beg);
			}
		}

		// read fine mask image data
		for (i = 0; i < ctrl->NUM_ROWS; i++) {
			mfdataFStream.read((char*)sensor_p[p]->fimage.mask[i], ctrl->NUM_COLS);
		}

		// read fine image data
		inval_fpc = 0;
		for (b = 0; b < ctrl->NUM_BANDS; b++) {
			for (i = 0; i < ctrl->NUM_ROWS; i++) {
				fdataFStream.read((char*)sensor_p[p]->fimage.data[b][i], ctrl->NUM_COLS * 2);
				// check invalid pixels
				for (j = 0; j < ctrl->NUM_COLS; j++) {
					if ((sensor_p[p]->fimage.data[b][i][j] < 0) ||								// negative reflectance -> invalid
						(sensor_p[p]->fimage.data[b][i][j] > sensor_p[p]->fimage.scale) ||		// reflectance > 1, i.e. > 1*scale_factor -> invalid
						(sensor_p[p]->fimage.data[b][i][j] == sensor_p[p]->fimage.fillv) ||		// invalid fillvalue -> invalid
						(sensor_p[p]->fimage.mask[i][j] > 0 )) {								// masked out -> invalid
							sensor_p[p]->fimage.data[b][i][j] = ctrl->INVALID;
							inval_fpc++;
					}
				}
			}
		}

		// read coarse image data
		inval_cpc = 0;
		if (!ctrl->CO_REGISTER.compare("NO")) {
			for (b = 0; b < ctrl->NUM_BANDS; b++) {
				for (i = 0; i < ctrl->NUM_ROWS; i++) {
					cdataFStream.read((char*)sensor_p[p]->cimage.data[b][i], ctrl->NUM_COLS * 2);
					// check invalid pixels
					for (j = 0; j < ctrl->NUM_COLS; j++) {
						if ((sensor_p[p]->cimage.data[b][i][j] < 0) ||								// negative reflectance -> invalid
							(sensor_p[p]->cimage.data[b][i][j] > sensor_p[p]->cimage.scale) || 		// reflectance > 1, i.e. > 1*scale_factor -> invalid
							(sensor_p[p]->cimage.data[b][i][j] == sensor_p[p]->cimage.fillv)) {		// invalid fillvalue -> invalid
							//(sensor_p[p]->cimage.mask[i][j] > 0)) {								// masked out -> invalid
							sensor_p[p]->cimage.data[b][i][j] = ctrl->INVALID;
							inval_cpc++;
						}
					}
				}
			}
		}
		else {
			// read coarse image data to be co-registered
			for (b = 0; b < ctrl->NUM_BANDS; b++) {
				for (i = 0; i < ctrl->COARSE_ROWS; i++) {
					cdataFStream.read((char*)cdata[b][i], ctrl->COARSE_COLS * 2);
					// check invalid pixels
					for (j = 0; j < ctrl->COARSE_COLS; j++) {
						if ((cdata[b][i][j] < 0) ||									// negative reflectance -> invalid
							(cdata[b][i][j] > sensor_p[p]->cimage.scale) || 		// reflectance > 1, i.e. > 1*scale_factor -> invalid
							(cdata[b][i][j] == sensor_p[p]->cimage.fillv)) {		// invalid fillvalue -> invalid
							//(sensor_p[p]->cimage.mask[i][j] > 0)) {				// masked out -> invalid
							cdata[b][i][j] = ctrl->INVALID;
							inval_cpc++;
						}
					}
				}
			}
			// perform co-registration and save the result
			string out_fname = sensor_p[p]->cimage.fname;
			size_t ext = out_fname.find(".dat");
			if (ext != string::npos)
			{
				//out_fname.erase(ext, 4);
				out_fname.insert(ext-12, "_matched");
			}
			else {
				//out_fname.append("_matched.dat");
				cout << "problem with co-registration output file name!" << endl;
			}
			co_register(sensor_p[p]->fimage.data, cdata, sensor_p[p]->cimage.data, ctrl, out_fname);
			log_file << "MODIS Co-register row = " << ctrl->CO_REGISTER_ROW << endl;
			log_file << "MODIS Co-register col = " << ctrl->CO_REGISTER_COL << endl;
		}

		fdataFStream.close();
		cdataFStream.close();
		mfdataFStream.close();

		cout << "Invalid pixels of all bands in fine image: " << inval_fpc << endl;
		cout << "Invalid pixels of all bands in coarse image: " << inval_cpc << endl;
		log_file << "Invalid pixels of all bands in fine image: " << inval_fpc << endl;
		log_file << "Invalid pixels of all bands in coarse image: " << inval_cpc << endl;
	}

	// free memory
	for (int b = 0; b < ctrl->NUM_BANDS; b++) {
		for (int i = 0; i < ctrl->COARSE_ROWS; i++) {
			delete[] cdata[b][i];
		}
		delete[] cdata[b];
	}
	delete[] cdata;

	return 0;
}

int PSRFM::read_pimage_data(short ***pdata, string fname, CONTROL_PARAMETER *ctrl) {

	int i;
	int j;
	int b;
	int bandSize;

	ifstream pdataFStream;  // prediction data file

	int num_bands;

	short ***cdata;

	// allocate memory for co-registration
	cdata = new short**[ctrl->NUM_BANDS];
	for (int b = 0; b < ctrl->NUM_BANDS; b++) {
		cdata[b] = new short*[ctrl->COARSE_ROWS];
		for (int r = 0; r < ctrl->COARSE_ROWS; r++) {
			cdata[b][r] = new short[ctrl->COARSE_COLS];
		}
	}

	if (!ctrl->CO_REGISTER.compare("NO")) {
		bandSize = ctrl->NUM_ROWS*ctrl->NUM_COLS;
	}
	else {
		bandSize = ctrl->COARSE_ROWS*ctrl->COARSE_COLS;
	}

	pdataFStream.open(fname.c_str(), ios::in | ios::binary | ios::ate);

	if (!pdataFStream.is_open()) {
		cout << "Cannot open input image file " << fname.c_str() << endl;
		return 1;
	}
	else { // read data
		num_bands = (int)pdataFStream.tellg() * sizeof(char) / (2 * bandSize);
		if (num_bands != ctrl->NUM_BANDS) {
			cout << "File not in correct dimension <" << fname.c_str() << ">" << endl;
		}
		else {
			// move to start of input File
			pdataFStream.seekg(0, ios_base::beg);
		}
	}

	if (!ctrl->CO_REGISTER.compare("NO")) {
		for (b = 0; b < ctrl->NUM_BANDS; b++) {
			for (i = 0; i < ctrl->NUM_ROWS; i++) {
				pdataFStream.read((char*)pdata[b][i], ctrl->NUM_COLS * 2);
			}
		}
	}
	else {
		for (b = 0; b < ctrl->NUM_BANDS; b++) {
			for (i = 0; i < ctrl->COARSE_ROWS; i++) {
				pdataFStream.read((char*)cdata[b][i], ctrl->COARSE_COLS * 2);
			}
		}
		// copy the co-registered part of the coarse image
		for (b = 0; b < ctrl->NUM_BANDS; b++) {
			for (i = 0; i < ctrl->NUM_ROWS; i++) {
				for (j = 0; j < ctrl->NUM_COLS; j++) {
					pdata[b][i][j] = cdata[b][i + ctrl->CO_REGISTER_ROW][j + ctrl->CO_REGISTER_COL];
				}
			}
		}
	}

	pdataFStream.close();

	// free memory
	for (int b = 0; b < ctrl->NUM_BANDS; b++) {
		for (int i = 0; i < ctrl->COARSE_ROWS; i++) {
			delete[] cdata[b][i];
		}
		delete[] cdata[b];
	}
	delete[] cdata;

	return 0;
}


int PSRFM::get_inputs(string ifname, SENSOR_PAIR *sensor_p[], string *pdata_fname, CONTROL_PARAMETER *ctrl, IMGAG_DATE_INFO *image_date_info, int PAIRS) {

	ifstream input;

	int i;
	string line;
	string word;
	int count = 0;

	vector<string> words;    // Create vector to hold our words

	int wordsCount;
	int rowss;
	int colss;

	input.open(ifname.c_str(), ios::in);

	if (!input.is_open()) {
		cout << "Error: Failed to open file <<" << ifname.c_str() << endl;
		return 1;
	}

	while (getline(input, line))
	{
		cout << line << '\n';

		if (!line.find("PSRFM_PARAMETER_START")) {
			cout << line.c_str() << endl;
		}
		else if (!line.find("#")) {
			cout << "Ignore... " << endl;
		}
		else {
			istringstream iss(line);
			words.clear();

			while (iss >> word) {
				words.push_back(word);
			}
			wordsCount = static_cast<int>(words.size());

			if (words[0] == "IN_PAIR_COARSE_FNAME") {
				cout << "coarse images..." << endl;
				for (i = 2; i < wordsCount; i++) {
					sensor_p[i - 2]->cimage.fname = words[i];

					get_image_date_info(words[i], &(image_date_info->date_num[i - 2]), &(image_date_info->c_date[i - 2]));


					cout << "..." << words[i].c_str() << endl;
				}
			}
			else if (words[0] == "IN_PAIR_FINE_FNAME") {
				cout << "fine images..." << endl;
				for (i = 2; i < wordsCount; i++) {
					sensor_p[i - 2]->fimage.fname = words[i];

					get_image_date_info(words[i], &(image_date_info->date_num[i - 2]), &(image_date_info->f_date[i - 2]));

					cout << "..." << words[i].c_str() << endl;
				}
			}
			else if (words[0] == "IN_PAIR_FINE_MASK_FNAME") {
				cout << "fine mask images..." << endl;
				for (i = 2; i < wordsCount; i++) {
					sensor_p[i - 2]->fimage.mfname = words[i];
					cout << "..." << words[i].c_str() << endl;
				}
			}
			else if (words[0] == "IN_PDAY_COARSE_FNAME") {
				cout << "Coarse image for prediction..." << endl;
				ctrl->NUM_PREDICTIONS = wordsCount - 2;

				for (i = 2; i < wordsCount; i++) {

					pdata_fname[i - 2] = words[i];
					int p = i - 2;

					get_image_date_info(words[i], &(image_date_info->pdt_date_num[p]), &(image_date_info->pdt_date[p]));

					cout << "..." << words[i].c_str() << endl;
				}
			}
			else if (words[0] == "OUT_PREDICTION_DIR") {
				cout << "OUTPUT_PREDICTION_DIR..." << endl;
				ctrl->OUTPUT_DIR = words[2];
			}
			else if (words[0] == "OUT_TEMP_DIR") {
				cout << "OUTPUT_TEMP_DIR..." << endl;
				ctrl->TEMP_DIR = words[2];
			}
			else if (words[0] == "OUT_ENVI_HDR") {
				cout << "OUTPUT_ENVI_HEDAER..." << endl;
				ctrl->ENVI_HDR = words[2];
			}
			else if (words[0] == "CO_REGISTER") {
				cout << "CO_REGISTER..." << endl;
				ctrl->CO_REGISTER = words[2];
			}
			else if (words[0] == "NROWS") {
				cout << "Number of rwos ...";
				rowss = stoi(words[2]);
				ctrl->NUM_ROWS = rowss;
				cout << rowss << endl;
			}
			else if (words[0] == "NCOLS") {
				cout << "Number of cols ...";
				colss = stoi(words[2]);
				ctrl->NUM_COLS = colss;
				cout << colss << endl;
			}
			else if (words[0] == "COARSE_ROWS") {
				cout << "Row number of coarse image ...";
				rowss = stoi(words[2]);
				ctrl->COARSE_ROWS = rowss;
				cout << rowss << endl;
			}
			else if (words[0] == "COARSE_COLS") {
				cout << "Colum number of coarse image ...";
				colss = stoi(words[2]);
				ctrl->COARSE_COLS = colss;
				cout << colss << endl;
			}
			else if (words[0] == "NBANDS") {
				cout << "Number of bands ...";
				int bands = stoi(words[2]);
				ctrl->NUM_BANDS = bands;
				cout << bands << endl;
			}
			else if (words[0] == "RESOLUTION") {
				cout << "Resolution ...";
				int resolution = stoi(words[2]);
				ctrl->RESOLUTION = resolution;
				cout << resolution << endl;
			}
			else if (words[0] == "BLOCK_SIZE") {
				cout << "Block size ...";
				int block_size = stoi(words[2]);
				ctrl->BLOCK_SIZE = block_size;
				cout << block_size << endl;
			}
			else if (words[0] == "SCALE_FACTOR") {
				cout << "SCALE_FACTOR ...";
				int scale_factor = stoi(words[2]);
				cout << scale_factor << endl;
			}
			else if (words[0] == "FINE_IMAGE_DATA_RANGE") {
				cout << "FINE_IMAGE_DATA_RANGE ...";
				int fine_image_data_min = stoi(words[2]);
				int fine_image_data_max = stoi(words[3]);
				for (i = 0; i < PAIRS; i++) {
					sensor_p[i]->fimage.min = fine_image_data_min;
					sensor_p[i]->fimage.max = fine_image_data_max;
				}
				cout << "MIN = " << fine_image_data_min << "  MAX = " << fine_image_data_max << endl;
			}
			else if (words[0] == "COARSE_IMAGE_DATA_RANGE") {
				cout << "COARSE_IMAGE_DATA_RANGE ...";
				int coarse_image_data_min = stoi(words[2]);
				int coarse_image_data_max = stoi(words[3]);
				for (int i = 0; i < PAIRS; i++) {
					sensor_p[i]->cimage.min = coarse_image_data_min;
					sensor_p[i]->cimage.max = coarse_image_data_max;
				}
				cout << "MIN = " << coarse_image_data_min << "  MAX = " << coarse_image_data_max << endl;
			}
			else if (words[0] == "FINE_IMAGE_FILLV") {
				cout << "FINE_IMAGE_FILLV ...";
				short fine_image_fill = short(stoi(words[2]));
				for (i = 0; i < PAIRS; i++) {
					sensor_p[i]->fimage.fillv = fine_image_fill;
				}
				cout << fine_image_fill << endl;
			}
			else if (words[0] == "COARSE_IMAGE_FILLV") {
				cout << "COARSE_IMAGE_FILLV ...";
				short coarse_image_fill = short(stoi(words[2]));
				for (i = 0; i < PAIRS; i++) {
					sensor_p[i]->fimage.fillv = coarse_image_fill;
				}
				cout << coarse_image_fill << endl;
			}

			else if (words[0] == "CLUSTER_RANGE") {
				cout << "CLUSTER_RANGE ...";
				ctrl->CLUSTER_RANGE[0] = stoi(words[2]);
				ctrl->CLUSTER_RANGE[1] = stoi(words[3]);
				cout << "MIN = " << ctrl->CLUSTER_RANGE[0] << "  MAX = " << ctrl->CLUSTER_RANGE[1] << endl;
			}
			else if (words[0] == "CLUSTER_METHOD") {
				cout << "CLUSTER_Method ...";
				ctrl->CLUSTER_METHOD = words[2];
				cout << ctrl->CLUSTER_METHOD << endl;
			}
			else if (words[0] == "CLUSTER_DATA") {
			cout << "CLUSTER_input data ...";
			ctrl->CLUSTER_DATA = words[2];
			cout << ctrl->CLUSTER_DATA << endl;
			}
			else if (words[0] == "FINE_IMAGE_UNCERTAINTY") {
				cout << "Fine uncertainty ...";
				ctrl->f_error = stof(words[2]);
				cout << ctrl->f_error << endl;
			}
			else if (words[0] == "COARSE_IMAGE_UNCERTAINTY") {
				cout << "Coarse uncertainty ...";
				ctrl->c_error = stof(words[2]);
				cout << ctrl->c_error << endl;
			}
			else if (words[0] == "CLUSTER_OPTIMAL") {
				cout << "CLUSTER_OPTIMAL...";
				ctrl->CLUSTER_OPTIMAL = words[2];
				cout << ctrl->CLUSTER_OPTIMAL << endl;
			}
			else if (words[0] == "RC_METHOD") {
			cout << "RC_METHOD...";
			ctrl->RC_METHOD = words[2];
			cout << ctrl->RC_METHOD << endl;
			}
			else if (words[0] == "MERGE_METHOD") {
				cout << "MERGE_METHOD ...";
				ctrl->MERGE_METHOD = words[2];
				cout << ctrl->MERGE_METHOD << endl;
			}
		}
	}

	input.close();

	// allocate data
	for (int p = 0; p < PAIRS; p++) {
		sensor_p[p]->fimage.data = new short**[ctrl->NUM_BANDS];
		sensor_p[p]->cimage.data = new short**[ctrl->NUM_BANDS];
		for (int b = 0; b < ctrl->NUM_BANDS; b++) {
			sensor_p[p]->fimage.data[b] = new short*[ctrl->NUM_ROWS];
			sensor_p[p]->cimage.data[b] = new short*[ctrl->NUM_ROWS];
			sensor_p[p]->fimage.mask = new unsigned char *[ctrl->NUM_ROWS];
			for (int r = 0; r < ctrl->NUM_ROWS; r++) {
				sensor_p[p]->fimage.data[b][r] = new short[ctrl->NUM_COLS];
				sensor_p[p]->cimage.data[b][r] = new short[ctrl->NUM_COLS];
				sensor_p[p]->fimage.mask[r] = new unsigned char[ctrl->NUM_COLS];
			}
		}
	}

	return 0;
}


void PSRFM::block_mean(short **f_idata, float **cf_idata, short **f_pdata, float **cf_pdata, int rows, int cols, int block_rows, int block_cols, int inval) {
	
	int i;
	int j;
	int m;
	int n;

	int ii = 0;
	int jj = 0;
	int inval_count = 0;

	int noPixs = block_rows * block_cols;

	for (i = 0; i < rows; i = i + block_rows) {
		for (j = 0; j < cols; j = j + block_cols) {

			cf_idata[ii][jj] = 0.0;
			cf_pdata[ii][jj] = 0.0;
			inval_count = 0;

			for (m = i; m < i + block_rows; m++) {
				for (n = j; n < j + block_cols; n++) {
					if (f_idata[m][n] == inval || f_pdata[m][n] == inval) {
						inval_count++;
					} 
					else {
						cf_idata[ii][jj] += f_idata[m][n];
						cf_pdata[ii][jj] += f_pdata[m][n];
					}
				}
			}

			if (inval_count == 0) {
				cf_idata[ii][jj] /= noPixs;
				cf_pdata[ii][jj] /= noPixs;
			}
			else {
				cf_idata[ii][jj] = float(inval);
				cf_pdata[ii][jj] = float(inval);
			}

			jj++;
		}

		ii++;
		jj = 0;
	}
}

int PSRFM::build_equations(short **f_idata, float **cf_idata, short **f_pdata, float **cf_pdata, unsigned char **clusters, float **A, float *l, float dt, int k, CONTROL_PARAMETER *ctrl) {

	int mk;

	int i, j, m, n, ik;

	int ii = 0;
	int jj = 0;
	int inval_count = 0;

	int noPixs = ctrl->BLOCK_SIZE * ctrl->BLOCK_SIZE;

	int *cluster_count;

	cluster_count = new int[k];

	mk = 0;

	for (i = 0; i < ctrl->NUM_ROWS; i = i + ctrl->BLOCK_SIZE) {
		for (j = 0; j < ctrl->NUM_COLS; j = j + ctrl->BLOCK_SIZE) {

			cf_idata[ii][jj] = 0.0;
			cf_pdata[ii][jj] = 0.0;
			inval_count = 0;

			for (ik = 0; ik < k; ik++) {
				cluster_count[ik] = 0;
			}

			for (m = i; m < i + ctrl->BLOCK_SIZE; m++) {
				for (n = j; n < j + ctrl->BLOCK_SIZE; n++) {
					if (f_idata[m][n] == ctrl->INVALID || f_pdata[m][n] == ctrl->INVALID) {
						inval_count++;
					}
					else {
						cf_idata[ii][jj] += f_idata[m][n];
						cf_pdata[ii][jj] += f_pdata[m][n];
						cluster_count[clusters[m][n]]++;
					}
				}
			}

			if (inval_count == 0) {
				cf_idata[ii][jj] /= noPixs;
				cf_pdata[ii][jj] /= noPixs;
				for (ik = 0; ik < k; ik++) {
					A[mk][ik] = float(cluster_count[ik]) / noPixs;					
				}
				l[mk] = (cf_pdata[ii][jj] - cf_idata[ii][jj]) / dt;
				mk++;
			}
			else {
				cf_idata[ii][jj] = float(ctrl->INVALID);
				cf_pdata[ii][jj] = float(ctrl->INVALID);
				for (ik = 0; ik < k; ik++) {
					A[mk][ik] = 0.0;
				}
				l[mk] = 0.0;
				mk++;
			}

			jj++;
		}

		ii++;
		jj = 0;
	}

	return mk;
}



// block means of input and prediction coarse resolution images and their changes for each coarse pixel
void PSRFM::block_mean_and_cc_change(float **oidata, short **idata, float **opdata, short **pdata, float *cc_change, int rows, int cols, int block_rows, int block_cols) {
	int i;
	int j;
	int m;
	int n;
	int k;

	int ii = 0;
	int jj = 0;
	k = 0;

	float numPixs = float(block_rows * block_cols);

	for (i = 0; i < rows; i = i + block_rows) {
		for (j = 0; j < cols; j = j + block_cols) {
			oidata[ii][jj] = 0.0;
			opdata[ii][jj] = 0.0;

			for (m = i; m < i + block_rows; m++) {
				for (n = j; n < j + block_cols; n++) {
					oidata[ii][jj] += idata[m][n];
					opdata[ii][jj] += pdata[m][n];
				}
			}
			oidata[ii][jj] /= numPixs;
			opdata[ii][jj] /= numPixs;
			cc_change[k] = opdata[ii][jj] - oidata[ii][jj];
			jj++;
			k++;
		}
		ii++;
		jj = 0;
	}
}

// blend_dir	forward or backward
// snesor_p		sensor data pairs
// pimage		prediction images
// pidx			prediction index 0, 1, ..., n predictions
int PSRFM::psrfm_blending(int blend_dir, int pidx, SENSOR_PAIR *sensor_p[], short ***pimage, CONTROL_PARAMETER *ctrl, IMGAG_DATE_INFO *image_date_info, short ***blend_output, short ***blend_uncertainty, string fine_sensor_name, ofstream& log_file) {

//	clock_t begin;
//	clock_t end;
//	double elapsed_secs;
	
	int k;							// for cluster loop
	int b;							// for band loop
	long inval_count = 0;			

	short ***in_cluster_data;
	int in_cluster_data_band;

	unsigned char **iClusters;
	unsigned char **saved_iClusters;
								
	short ***saved_blend_image;
	unsigned short ***saved_blend_image_uncertainty;
	

	ofstream prd_output;			// prediction output stream
	string prd_output_fname;		// prediction output file name

	ofstream unctn_output;			// uncertainty output stream
	string unctn_output_fname;		// uncertainty output file name

	ofstream cluster_output;		// cluster output stream
	string cluster_output_fname;	// cluster output file name

	int ORG_C_ROWS;					// original number of rows of coarse image 
	int ORG_C_COLS;					// original number of columns of coarse image 

	ORG_C_ROWS = ctrl->NUM_ROWS / ctrl->BLOCK_SIZE;
	ORG_C_COLS = ctrl->NUM_COLS / ctrl->BLOCK_SIZE;

	// for quality assessment
	double CC = 0.0;					// correlation at coarse scale of fine difference and coarse difference
	double CF = 0.0;					// correlation at fine scale of fine difference and coarse difference
	double RR = 0.0;
	double Sigma = 0.0;

	double sum_CC = 0.0;				// correlation at coarse scale of fine difference and coarse difference
	double sum_CF = 0.0;				// correlation at fine scale of fine difference and coarse difference
	double sum_RR = 0.0;
	double sum_Sigma = 0.0;

	double CC2 = 0.0;					// correlation at coarse scale of fine difference and coarse difference
	double CF2 = 0.0;					// correlation at fine scale of fine difference and coarse difference
	double RR2 = 0.0;
	double Sigma2 = 0.0;

	double sum_CC2 = 0.0;				// correlation at coarse scale of fine difference and coarse difference
	double sum_CF2 = 0.0;				// correlation at fine scale of fine difference and coarse difference
	double sum_RR2 = 0.0;
	double sum_Sigma2 = 0.0;

	double cluster_optimal_value;
	double cluster_optimal_value_old;


	int max_k = 0;

	float **fc_idata;
	float **fc_pdata;

	int tmp_count = 1;
	double RS2;
	double RS2_sum = 0.0;
	double RS2_ave = 0.0;

	float **dcc;				// different coarse images at coarse resolution
	float **dfc;				// different fine image at coarse resolution

	float **rs;					// residuals
	float **rsc;				// residual corrections
	short **rsc_blend;			// residual corrections adjusted blend image

	float **dff;				// different fine at fine resolution
	float **dcf;				// different coarse at fine resolution

	dff = new float *[ctrl->NUM_ROWS];
	dcf = new float *[ctrl->NUM_ROWS];
	rsc = new float *[ctrl->NUM_ROWS];
	rsc_blend = new short *[ctrl->NUM_ROWS];
	for (int i = 0; i < ctrl->NUM_ROWS; i++) {
		dff[i] = new float[ctrl->NUM_COLS];
		dcf[i] = new float[ctrl->NUM_COLS];
		rsc[i] = new float[ctrl->NUM_COLS];
		rsc_blend[i] = new short[ctrl->NUM_COLS];
	}

	fc_idata = new float *[ORG_C_ROWS];
	fc_pdata = new float *[ORG_C_ROWS];
	dcc = new float *[ORG_C_ROWS];
	dfc = new float *[ORG_C_ROWS];
	rs = new float *[ORG_C_ROWS];
	for (int i = 0; i < ORG_C_ROWS; i++) {
		fc_idata[i] = new float[ORG_C_COLS];
		fc_pdata[i] = new float[ORG_C_COLS];
		dcc[i] = new float[ORG_C_COLS];
		dfc[i] = new float[ORG_C_COLS];
		rs[i] = new float[ORG_C_COLS];
	}

	saved_blend_image = new short **[ctrl->NUM_BANDS];
	saved_blend_image_uncertainty = new unsigned short **[ctrl->NUM_BANDS];
	for (b = 0; b < ctrl->NUM_BANDS; b++) {
		saved_blend_image[b] = new short *[ctrl->NUM_ROWS];
		saved_blend_image_uncertainty[b] = new unsigned short *[ctrl->NUM_ROWS];
		for (int i = 0; i < ctrl->NUM_ROWS; i++) {
			saved_blend_image[b][i] = new short [ctrl->NUM_COLS];
			saved_blend_image_uncertainty[b][i] = new unsigned short[ctrl->NUM_COLS];
		}
	}

	if (!ctrl->CLUSTER_DATA.compare("fine+coarse")) {
		in_cluster_data = new short **[ctrl->NUM_BANDS * 2];
		for (b = 0; b < ctrl->NUM_BANDS * 2; b++) {
			in_cluster_data[b] = new short *[ctrl->NUM_ROWS];
			for (int i = 0; i < ctrl->NUM_ROWS; i++) {
				in_cluster_data[b][i] = new short[ctrl->NUM_COLS];
			}
		}
	}
	else {
		in_cluster_data = new short **[ctrl->NUM_BANDS];
		for (b = 0; b < ctrl->NUM_BANDS; b++) {
			in_cluster_data[b] = new short *[ctrl->NUM_ROWS];
			for (int i = 0; i < ctrl->NUM_ROWS; i++) {
				in_cluster_data[b][i] = new short[ctrl->NUM_COLS];
			}
		}
	}

	saved_iClusters = new unsigned char *[ctrl->NUM_ROWS];
	for (int i = 0; i < ctrl->NUM_ROWS; i++) {
		saved_iClusters[i] = new unsigned char[ctrl->NUM_COLS];
	}

	int mk;				// total number of coarse resolution image pixels;

	int BLOCK_ROWS = ctrl->BLOCK_SIZE;
	int BLOCK_COLS = ctrl->BLOCK_SIZE;

	float dt;			// temporal interval of input and prediction images
	int nbThreads;

		
	mk = ORG_C_ROWS * ORG_C_COLS;

	
	short **blended_image;
	blended_image = new short *[ctrl->NUM_ROWS];
	for (int i = 0; i < ctrl->NUM_ROWS; i++) {
		blended_image[i] = new short[ctrl->NUM_COLS];
	}


	float **A;				// cluster fraction of a block (a coarse resolution image pixel)
	float *cc_refl_change;


	cc_refl_change = new float[ORG_C_ROWS * ORG_C_COLS];

	
	float **cc_idata = new float *[ORG_C_ROWS];		// coarse image in coarse resoluton
	float **cc_pdata = new float *[ORG_C_ROWS];		// coarse prediction image in coarse resolution
	for (int i = 0; i < ORG_C_ROWS; i++) {
		cc_idata[i] = new float[ORG_C_COLS];
		cc_pdata[i] = new float[ORG_C_COLS];
	}


	// for OpenMP
	// get number of cores and use all but minus 1
	int numCores = omp_get_num_procs();
	numCores--;

	Eigen::initParallel();
	nbThreads = Eigen::nbThreads() - 1;
	setNbThreads(nbThreads);


	// output file name is sensor name of fine image + like + prediction date + fwd (bwd) + .dat
	prd_output_fname = ctrl->TEMP_DIR;
	prd_output_fname.append("\\");
	prd_output_fname.append(fine_sensor_name);
	prd_output_fname.append("_like_");
	prd_output_fname.append(image_date_info->pdt_date[pidx]);

	unctn_output_fname = prd_output_fname;
	
	int im;			// im = 0 forward. im = 1 backward

	if (blend_dir == FORWARD) {
		im = 0;
		prd_output_fname.append("_FWD.dat");
		unctn_output_fname.append("_FWD_Q.dat");
	}
	else if (blend_dir == BACKWARD) {
		im = 1;
		prd_output_fname.append("_BWD.dat");
		unctn_output_fname.append("_BWD_Q.dat");
	}

	iClusters = new unsigned char *[ctrl->NUM_ROWS];
	
	for (int i = 0; i < ctrl->NUM_ROWS; i++) {
		iClusters[i] = new unsigned char[ctrl->NUM_COLS];
		
	}

	prd_output.open(prd_output_fname.c_str(), ios::out | ios::binary);
	if (!prd_output) {
		cout << "Cannot open file for writing: <" << prd_output_fname.c_str() << ">" << endl;
	}

	unctn_output.open(unctn_output_fname.c_str(), ios::out | ios::binary);
	if (!unctn_output) {
		cout << "Cannot open file for writing: <" << unctn_output_fname.c_str() << ">" << endl;
	}

	k = ctrl->CLUSTER_RANGE[0];		// initial number of clusters

	
	MatrixXf EigenA = MatrixXf::Zero(mk, k);			// A Ax = b
	VectorXf Eigenlv = VectorXf::Zero(mk);				// b
	MatrixXf EigenA2 = MatrixXf::Zero(mk, k);			// weighted A
	VectorXf Eigenlv2 = VectorXf::Zero(mk);				// weighted lv

	//VectorXf EigenB = VectorXf::Zero(k);				// A'Pb

	VectorXf Eigenfv = VectorXf::Zero(k);				// x
	MatrixXf EigenPv = MatrixXf::Zero(mk, mk);			// weights

	MatrixXf EigenNv = MatrixXf::Zero(k, k);			// k cluster number
	MatrixXf EigenQv = MatrixXf::Zero(mk, mk);
	MatrixXf EigenQvv = MatrixXf::Zero(mk, mk);
	//MatrixXf EigenRv0 = MatrixXf::Zero(mk, mk);
	MatrixXf EigenRv = MatrixXf::Zero(mk, mk);
	MatrixXf EigenRvDiag = MatrixXf::Zero(mk, mk);

	VectorXf Eigenvv = VectorXf::Zero(mk);

	float sv2;
	float rv;


	// time interval between the input and the prediction date
	dt = float(abs(image_date_info->date_num[im] - image_date_info->pdt_date_num[pidx]));

	cluster_optimal_value_old = 1000000.0;				// cluster optimal value holder

	// input data for clustering 
	if (!ctrl->CLUSTER_DATA.compare("fine")) {			// fine image only
		// copy the fine image data
		for (b = 0; b < ctrl->NUM_BANDS; b++) {
			for (int i = 0; i < ctrl->NUM_ROWS; i++) {
				for (int j = 0; j < ctrl->NUM_COLS; j++) {
					in_cluster_data[b][i][j] = short(sensor_p[im]->fimage.data[b][i][j]);
				}

			}
		}
		in_cluster_data_band = ctrl->NUM_BANDS;
	}
	else if (!ctrl->CLUSTER_DATA.compare("fine+coarse")) {	// fine + coarse image	
		// copy the fine image data
		for (b = 0; b < ctrl->NUM_BANDS; b++) {
			for (int i = 0; i < ctrl->NUM_ROWS; i++) {
				for (int j = 0; j < ctrl->NUM_COLS; j++) {
					in_cluster_data[b][i][j] = short(sensor_p[im]->fimage.data[b][i][j]);
				}

			}
		}
		// append the coarse image data at prediction date
		for (b = 0; b < ctrl->NUM_BANDS; b++) {
			for (int i = 0; i < ctrl->NUM_ROWS; i++) {
				for (int j = 0; j < ctrl->NUM_COLS; j++) {
					in_cluster_data[b + ctrl->NUM_BANDS][i][j] = short(pimage[b][i][j]);
				}

			}
		}
		in_cluster_data_band = ctrl->NUM_BANDS*2;
	}
	else if (!ctrl->CLUSTER_DATA.compare("ratio")) {
		// derive the change ratios
		for (b = 0; b < ctrl->NUM_BANDS; b++) {
			for (int i = 0; i < ctrl->NUM_ROWS; i++) {
				for (int j = 0; j < ctrl->NUM_COLS; j++) {
					// deal with invalid pixels
					if (pimage[b][i][j] == ctrl->INVALID ||
						sensor_p[im]->cimage.data[b][i][j] == ctrl->INVALID ||
						sensor_p[im]->cimage.data[b][i][j] == 0) {
						in_cluster_data[b][i][j] = ctrl->INVALID;
					}
					else {
						in_cluster_data[b][i][j] = short((pimage[b][i][j] - sensor_p[im]->cimage.data[b][i][j]) * sensor_p[im]->cimage.scale / sensor_p[im]->cimage.data[b][i][j]);
					}
				}

			}
		}
		in_cluster_data_band = ctrl->NUM_BANDS;
	}

	// set the first band pixel of the in_cluster_data to invalid if any one band is invalid
	// to save loop time for the clustering method KMEAN	
	if (!ctrl->CLUSTER_METHOD.compare("KMEAN")) {
		for (b = 0; b < in_cluster_data_band; b++) {
			for (int i = 0; i < ctrl->NUM_ROWS; i++) {
				for (int j = 0; j < ctrl->NUM_COLS; j++) {
					if (in_cluster_data[b][i][j] == ctrl->INVALID) {
						in_cluster_data[0][i][j] = ctrl->INVALID;
						inval_count++;
						break;
					}
				}
			}
		}
	}

	// for clusters 
	for (k = ctrl->CLUSTER_RANGE[0]; k <= ctrl->CLUSTER_RANGE[1]; k++) {
	
		cout << endl << "Cluster No: " << k << endl;
		log_file << endl << "Cluster No: " << k << endl;

		sum_CC = 0.0;
		sum_CF = 0.0;
		sum_RR = 0.0;
		sum_Sigma = 0.0;

		EigenA.resize(Eigen::NoChange, k);			
		EigenA2.resize(Eigen::NoChange, k);			
		Eigenfv.resize(k);							
		EigenNv.resize(k, k);
		//EigenB.resize(k);

		A = new float *[mk];
		for (int ik = 0; ik < mk; ik++) {
			A[ik] = new float[k];
		}

		// for each band
		bool cluster_done = false;
		for (b = 0; b < ctrl->NUM_BANDS; b++) {
		
			// clustering 
			if (!ctrl->CLUSTER_METHOD.compare("KMEAN") && !cluster_done) {		// kemen using all the bands for clustering
				if (!ctrl->CLUSTER_DATA.compare("fine+coarse")) {	// fine + coarse image
					KMeans(in_cluster_data, iClusters, k, ctrl->NUM_BANDS*2, ctrl->NUM_ROWS, ctrl->NUM_COLS, ctrl->ITERATIONS, ctrl->MAX_DIFF_THRESHOLD, ctrl->INVALID);
				}
				else {
					KMeans(in_cluster_data, iClusters, k, ctrl->NUM_BANDS, ctrl->NUM_ROWS, ctrl->NUM_COLS, ctrl->ITERATIONS, ctrl->MAX_DIFF_THRESHOLD, ctrl->INVALID);
				}
				
				cluster_done = true;

			}
			else if (!ctrl->CLUSTER_METHOD.compare("CRATIO")) {  // simple clustering using change ratio of single band
				CRatio(in_cluster_data, iClusters, k, b, ctrl->NUM_ROWS, ctrl->NUM_COLS, ctrl->INVALID);
			}
			
			/*
			// build observation equations, the following code was replaced by calling the function build_equations()
			block_mean(sensor_p[im]->cimage.data[b], cc_idata, pimage[b], cc_pdata, ctrl->NUM_ROWS, ctrl->NUM_COLS, BLOCK_ROWS, BLOCK_COLS, ctrl->INVALID);
			// get the cluster percetage within a coarse resolution pixel -- A is the output
			mk = area_fraction(iClusters, k, ctrl->NUM_ROWS, ctrl->NUM_COLS, ctrl->BLOCK_SIZE, A);
			// derive the change rate 
			int ii = 0;
			for (int i = 0; i < ORG_C_ROWS; i++) {
				for (int j = 0; j < ORG_C_COLS; j++) {
					// deal with invalid pixels
					if (cc_pdata[i][j] == ctrl->INVALID || cc_idata[i][j] == ctrl->INVALID) {
						cc_refl_change[ii] = 0;
					}
					else {
						cc_refl_change[ii] = (cc_pdata[i][j] - cc_idata[i][j]) / dt;
					}
					ii++;
				}
			}
			*/

			mk = build_equations(sensor_p[im]->cimage.data[b], cc_idata, pimage[b], cc_pdata, iClusters, A, cc_refl_change, dt, k, ctrl);

			// for weight based on uncertainty input for linear system equations
			double wu = sqrt(2.0)*ctrl->c_error / dt;
			wu = 1.0/(wu*wu);
			for (int i = 0; i < mk; i++) {
				EigenPv(i, i) = float(wu);
			}

			for (int i = 0; i < mk; i++) {
				Eigenlv(i) = cc_refl_change[i];		// change rate for Eigen operation
				for (int j = 0; j < k; j++) {
					EigenA(i, j) = A[i][j];
				}
			}

			// weighting with uncertainty
			//EigenA2.noalias() = EigenPv * EigenA;
			//Eigenlv2.noalias() = EigenPv * Eigenlv;
			EigenA2 = EigenPv * EigenA;
			Eigenlv2 = EigenPv * Eigenlv;

			// EIGEN method to solve linear systems of equations
			Eigenfv = EigenA2.bdcSvd(ComputeThinU | ComputeThinV).solve(Eigenlv2);
			EigenNv.noalias() = EigenA.transpose() * EigenPv * EigenA;
			EigenQv = EigenNv.inverse();

			//EigenNv.noalias() = EigenA.transpose() * EigenA2;
			//EigenNv = EigenA.transpose() * EigenA2;
			//EigenB.noalias() = EigenA.transpose() * Eigenlv2;
			//EigenB = EigenA.transpose() * Eigenlv2;
			//EigenQv = EigenNv.inverse();
			//EigenQv = EigenNv.completeOrthogonalDecomposition().pseudoInverse();			
			//Eigenfv = EigenQv * EigenB;

			//EigenRv.noalias() = EigenPv - EigenA2 * EigenQv * EigenA2.transpose();	// Rv
			EigenRv = EigenPv - EigenA2 * EigenQv * EigenA2.transpose();	// Rv
			Eigenvv = EigenA * Eigenfv - Eigenlv;			// vv
			EigenRvDiag = EigenRv.diagonal();
			rv = EigenRvDiag.sum();							// rv

			sv2 = Eigenvv.transpose()*EigenPv*Eigenvv;
			sv2 = sv2 / rv;									// sv2
			Sigma = sqrt(sv2);
			sum_Sigma += Sigma;

			EigenQvv = sv2 * EigenQv;						// Qvv


			// get the predicted image and its uncertainty
			for (int i = 0; i < ctrl->NUM_ROWS; i++) {
				for (int j = 0; j < ctrl->NUM_COLS; j++) {
					// deal with invalid pixels
					if (sensor_p[im]->fimage.data[b][i][j] == ctrl->INVALID || iClusters[i][j] > k-1) {	// cluster k+1 is used for invalid pixels
						blend_output[b][i][j] = short(ctrl->INVALID); // short(sensor_p[im]->fimage.fillv);
						blend_uncertainty[b][i][j] = short(ctrl->INVALID);
					}
					else {
						short predict = short(sensor_p[im]->fimage.data[b][i][j] + dt * Eigenfv(iClusters[i][j]));
						if (predict < 0 || predict > sensor_p[im]->fimage.scale) {
							blend_output[b][i][j] = sensor_p[im]->fimage.data[b][i][j];
							blend_uncertainty[b][i][j] = short(ctrl->f_error*ctrl->f_error);
						}
						else {
							blend_output[b][i][j] = predict;
							blend_uncertainty[b][i][j] = short(ctrl->f_error*ctrl->f_error + dt * dt*EigenQvv(iClusters[i][j], iClusters[i][j]));
						}
					}
				}				
			}

			// residual & correlation analysis

			block_mean(sensor_p[im]->fimage.data[b], fc_idata, blend_output[b], fc_pdata, ctrl->NUM_ROWS, ctrl->NUM_COLS, BLOCK_ROWS, BLOCK_COLS, ctrl->INVALID);

			tmp_count = 1;
			RS2_ave = 0.0;
			RS2_sum = 0.0;
			for (int ii = 0; ii < ORG_C_ROWS; ii++) {
				for (int jj = 0; jj < ORG_C_COLS; jj++) {					
					// coarse image difference (input and the input for prediction) in coarse resolution
					if (cc_idata[ii][jj] == ctrl->INVALID || cc_pdata[ii][jj] == ctrl->INVALID) {
						dcc[ii][jj] = float(ctrl->INVALID);		// deal with invalid pixels
					}
					else {
						dcc[ii][jj] = cc_idata[ii][jj] - cc_pdata[ii][jj];		
					}
					// fine image difference (input and the blended) in coarse resolution
					if (fc_idata[ii][jj] == ctrl->INVALID || fc_pdata[ii][jj] == ctrl->INVALID) {
						dfc[ii][jj] = float(ctrl->INVALID);
					}
					else {
						dfc[ii][jj] = fc_idata[ii][jj] - fc_pdata[ii][jj];		
					}
					// fine and coarse difference in coarse resolution -> residuales
					if (dcc[ii][jj] == ctrl->INVALID || dfc[ii][jj] == ctrl->INVALID) {
						rs[ii][jj] = float(ctrl->INVALID);
					}
					else {
						rs[ii][jj] = dcc[ii][jj] - dfc[ii][jj];					
						RS2 = rs[ii][jj] * rs[ii][jj];
						RS2_sum += RS2;
						RS2_ave += (RS2 - RS2_ave) / tmp_count;
						tmp_count++;
					}
				}
			}
			
			Sigma = sqrt(RS2_sum / tmp_count);
			sum_Sigma += Sigma;

			for (int i = 0; i < ctrl->NUM_ROWS; i++) {
				for (int j = 0; j < ctrl->NUM_COLS; j++) {
					// deal with invalid pixels
					if (blend_output[b][i][j] == ctrl->INVALID || sensor_p[im]->fimage.data[b][i][j] == ctrl->INVALID) {
						dff[i][j] = float(ctrl->INVALID);
					}
					else {
						dff[i][j] = float(blend_output[b][i][j] - sensor_p[im]->fimage.data[b][i][j]);
					}
					if (pimage[b][i][j] == ctrl->INVALID || sensor_p[im]->cimage.data[b][i][j] == ctrl->INVALID) {
						dcf[i][j] = float(ctrl->INVALID);
					}
					else {
						dcf[i][j] = float(pimage[b][i][j] - sensor_p[im]->cimage.data[b][i][j]);
					}	
				}
			}

			// correlation of fine image difference (input and the blended - from predicted fine image) (from predicted fine image) and MODIS image for prediction at coarse resolution
			CC = abs(correlation(dfc, dcc, ORG_C_ROWS, ORG_C_COLS, ctrl->INVALID));
			sum_CC += CC;

			// correlation of fine image and MODIS image for prediction
			CF = abs(correlation(dff, dcf, ctrl->NUM_ROWS, ctrl->NUM_COLS, ctrl->INVALID));
			sum_CF += CF;

			// residues
			RR = sqrt( RS2_ave );
			RR = sqrt(RS2_sum /tmp_count);
			sum_RR += RR;

			cout << "B = " << b << endl;
			log_file << "B = " << b << endl;

			// residual adjustments
			if (!ctrl->RC_METHOD.compare("biliear")) {
				
				// compute residual corrections
				residual_adjust(rs, rsc, ctrl->NUM_ROWS, ctrl->NUM_COLS, BLOCK_ROWS, BLOCK_COLS, ctrl->INVALID);
				
				// add residual corrections
				for (int i = 0; i < ctrl->NUM_ROWS; i++) {
					for (int j = 0; j < ctrl->NUM_COLS; j++) {
						float after_rsc = blend_output[b][i][j] - rsc[i][j];
						if (after_rsc > 0 || after_rsc < sensor_p[im]->fimage.scale) {
							rsc_blend[i][j] = short(after_rsc);
						}
					}
				}
				
				// re-evaluate
				block_mean(sensor_p[im]->fimage.data[b], fc_idata, rsc_blend, fc_pdata, ctrl->NUM_ROWS, ctrl->NUM_COLS, BLOCK_ROWS, BLOCK_COLS, ctrl->INVALID);

				tmp_count = 1;
				RS2_ave = 0.0;
				RS2_sum = 0.0;
				for (int ii = 0; ii < ORG_C_ROWS; ii++) {
					for (int jj = 0; jj < ORG_C_COLS; jj++) {
						// the following code is not necessary because dcc was calculated already.
						// coarse image difference (input and the input for prediction) in coarse resolution
						/*if (cc_idata[ii][jj] == ctrl->INVALID || cc_pdata[ii][jj] == ctrl->INVALID) {
							dcc[ii][jj] = float(ctrl->INVALID);		// deal with invalid pixels
						}
						else {
							dcc[ii][jj] = cc_idata[ii][jj] - cc_pdata[ii][jj];
						}*/
						// fine image difference (input and the blended) in coarse resolution
						if (fc_idata[ii][jj] == ctrl->INVALID || fc_pdata[ii][jj] == ctrl->INVALID) {
							dfc[ii][jj] = float(ctrl->INVALID);
						}
						else {
							dfc[ii][jj] = fc_idata[ii][jj] - fc_pdata[ii][jj];
						}
						// fine and coarse difference in coarse resolution -> residuales
						if (dcc[ii][jj] == ctrl->INVALID || dfc[ii][jj] == ctrl->INVALID) {
							rs[ii][jj] = float(ctrl->INVALID);
						}
						else {
							rs[ii][jj] = dcc[ii][jj] - dfc[ii][jj];
							RS2 = rs[ii][jj] * rs[ii][jj];
							RS2_sum += RS2;
							RS2_ave += (RS2 - RS2_ave) / tmp_count;
							tmp_count++;
						}
					}
				}
				Sigma2 = sqrt(RS2_sum / tmp_count);
				sum_Sigma2 += Sigma2;

				for (int i = 0; i < ctrl->NUM_ROWS; i++) {
					for (int j = 0; j < ctrl->NUM_COLS; j++) {
						// deal with invalid pixels
						if (rsc_blend[i][j] == ctrl->INVALID || sensor_p[im]->fimage.data[b][i][j] == ctrl->INVALID) {
							dff[i][j] = float(ctrl->INVALID);
						}
						else {
							dff[i][j] = float(rsc_blend[i][j] - sensor_p[im]->fimage.data[b][i][j]);
						}
						// the following code is not necessary because dcf was calculated already.
						/*if (pimage[b][i][j] == ctrl->INVALID || sensor_p[im]->cimage.data[b][i][j] == ctrl->INVALID) {
							dcf[i][j] = float(ctrl->INVALID);
						}
						else {
							dcf[i][j] = float(pimage[b][i][j] - sensor_p[im]->cimage.data[b][i][j]);
						}*/
					}
				}

				// correlation of fine image difference (input and the blended - from predicted fine image) (from predicted fine image) and MODIS image for prediction at coarse resolution
				CC2 = abs(correlation(dfc, dcc, ORG_C_ROWS, ORG_C_COLS, ctrl->INVALID));
				sum_CC2 += CC2;

				// correlation of fine image and MODIS image for prediction
				CF2 = abs(correlation(dff, dcf, ctrl->NUM_ROWS, ctrl->NUM_COLS, ctrl->INVALID));
				sum_CF2 += CF2;

				// residues
				RR2 = sqrt(RS2_ave);
				RR2 = sqrt(RS2_sum / tmp_count);
				sum_RR2 += RR2;

				cout << "S1 = " << Sigma << "; CC1 = " << CC << "; CF1 = " << CF << ", RR1 = " << RR << endl;
				log_file << "S1 = " << Sigma << "; CC1 = " << CC << "; CF1 = " << CF << ", RR1 = " << RR << endl;
				cout << "S2 = " << Sigma << "; CC2 = " << CC2 << "; CF2 = " << CF2 << ", RR2 = " << RR2 << endl;
				log_file << "S2 = " << Sigma << "; CC2 = " << CC2 << "; CF2 = " << CF2 << ", RR2 = " << RR2 << endl;

				// compare the sum of the residual squares before and after the residual corrections
				if (CF2 > CF && CC2 > CC && RR2 < RR) {
					CF = CF2;
					RR = RR2;
					CC = CC2;
					for (int i = 0; i < ctrl->NUM_ROWS; i++) {
						for (int j = 0; j < ctrl->NUM_COLS; j++) {
							blend_output[b][i][j] = rsc_blend[i][j];
						}
					}
				}
			}

			cout << "S = " << Sigma << "; CC = " << CC << "; CF = " << CF << ", RR = " << RR << endl;
			log_file << "S = " << Sigma << "; CC = " << CC << "; CF = " << CF << ", RR = " << RR << endl;

		}  // end of band loop

		if (!ctrl->CLUSTER_OPTIMAL.compare("CC")) {
			cluster_optimal_value = 1.0 - sum_CC / ctrl->NUM_BANDS;
		}
		else if (!ctrl->CLUSTER_OPTIMAL.compare("CF")) {
			cluster_optimal_value = 1.0 - sum_CF / ctrl->NUM_BANDS;
		}
		else if (!ctrl->CLUSTER_OPTIMAL.compare("RR")) {
			cluster_optimal_value = sum_RR;
		}

		if (cluster_optimal_value < cluster_optimal_value_old) {
			cluster_optimal_value_old = cluster_optimal_value;

			max_k = k;				// max value with k cluster

			for (int b = 0; b < ctrl->NUM_BANDS; b++) {
				for (int i = 0; i < ctrl->NUM_ROWS; i++) {
					for (int j = 0; j < ctrl->NUM_COLS; j++) {
						saved_blend_image[b][i][j] = blend_output[b][i][j];
						saved_blend_image_uncertainty[b][i][j] = (unsigned short) sqrt(blend_uncertainty[b][i][j]);
					}
				}
			}

			for (int i = 0; i < ctrl->NUM_ROWS; i++) {
				for (int j = 0; j < ctrl->NUM_COLS; j++) {
					saved_iClusters[i][j] = iClusters[i][j];
				}
			}
		}

		for (int i = 0; i < k; i++) {
			delete[] A[i];
		}
		delete[] A;

	} // end of k (clustering) loop

	if (!ctrl->CLUSTER_OPTIMAL.compare("CC")) {
		cout << endl << "Cluster number: " << max_k << " has the maximun average value: CC = " << sum_CC / ctrl->NUM_BANDS << endl;
		log_file << endl << "Cluster number: " << max_k << " has the maximun average value: CC = " << sum_CC / ctrl->NUM_BANDS << endl;
	}
	else if (!ctrl->CLUSTER_OPTIMAL.compare("CF")) {
		cout << endl << "Cluster number: " << max_k << " has the maximun average value: CF = " << sum_CF / ctrl->NUM_BANDS << endl;
		log_file << endl << "Cluster number: " << max_k << " has the maximun average value: CF = " << sum_CF / ctrl->NUM_BANDS << endl;
	}
	else if (!ctrl->CLUSTER_OPTIMAL.compare("RR")) {
		cout << endl << "Cluster number: " << max_k << " has the smallest average value: RR = " << sum_RR / ctrl->NUM_BANDS << endl;
		log_file << endl << "Cluster number: " << max_k << " has the smallest average value: RR = " << sum_RR / ctrl->NUM_BANDS << endl;
	}
	
	

	// output the optimized blended image
	for (int b = 0; b < ctrl->NUM_BANDS; b++) {
		for (int i = 0; i < ctrl->NUM_ROWS; i++) {
			prd_output.write(reinterpret_cast<char *>(saved_blend_image[b][i]), sizeof(short)*ctrl->NUM_COLS);
			unctn_output.write(reinterpret_cast<char *>(saved_blend_image_uncertainty[b][i]), sizeof(unsigned short)*ctrl->NUM_COLS);
		}
	}
	prd_output.close();
	unctn_output.close();

	cluster_output_fname = ctrl->TEMP_DIR; //    "D:\\cluster_";	// cluster output file name
	cluster_output_fname.append("\\");
	cluster_output_fname.append(fine_sensor_name);
	cluster_output_fname.append("_cluster_");
	cluster_output_fname.append(image_date_info->pdt_date[pidx]);
	if (blend_dir == FORWARD) {
		prd_output_fname.append(image_date_info->pdt_date[pidx]);
		cluster_output_fname.append("_FWD_");
	}
	else {
		cluster_output_fname.append("_BWD_");
	}

	cluster_output_fname.append(to_string(max_k));
	cluster_output_fname.append(".dat");

	cluster_output.open(cluster_output_fname.c_str(), ios::out | ios::binary);

	// output cluster of the best image
	for (int i = 0; i < ctrl->NUM_ROWS; i++) {
		cluster_output.write(reinterpret_cast<char *>(saved_iClusters[i]), sizeof(unsigned char)*ctrl->NUM_COLS);
	}
	cluster_output.close();

	// return the optimized blended images
	for (int b = 0; b < ctrl->NUM_BANDS; b++) {
		for (int i = 0; i < ctrl->NUM_ROWS; i++) {
			for (int j = 0; j < ctrl->NUM_COLS; j++) {
				blend_output[b][i][j] = saved_blend_image[b][i][j];
				blend_uncertainty[b][i][j] = saved_blend_image_uncertainty[b][i][j] * saved_blend_image_uncertainty[b][i][j];
			}
		}
	}

	// free memory
	if (!ctrl->CLUSTER_DATA.compare("fine+coarse")) {
		for (int b = 0; b < ctrl->NUM_BANDS*2; b++) {
			for (int i = 0; i < ctrl->NUM_ROWS; i++) {
				delete[] in_cluster_data[b][i];
			}
			delete[] in_cluster_data[b];
		}
	}
	else {
		for (int b = 0; b < ctrl->NUM_BANDS; b++) {
			for (int i = 0; i < ctrl->NUM_ROWS; i++) {
				delete[] in_cluster_data[b][i];
			}
			delete[] in_cluster_data[b];
		}
	}
	delete[] in_cluster_data;

	for (int b = 0; b < ctrl->NUM_BANDS; b++) {
		for (int i = 0; i < ctrl->NUM_ROWS; i++) {
			delete[] saved_blend_image[b][i];
			delete[] saved_blend_image_uncertainty[b][i];
		}
		delete[] saved_blend_image[b];
		delete[] saved_blend_image_uncertainty[b];
	}
	delete[] saved_blend_image;
	delete[] saved_blend_image_uncertainty;


	for (int i = 0; i < ctrl->NUM_ROWS; i++) {
		delete[] iClusters[i];
		delete[] dff[i];
		delete[] dcf[i];
		delete[] rsc[i];
		delete[] rsc_blend[i];
	}

	delete[] iClusters;
	delete[] dff;
	delete[] dcf;
	delete[] rsc;
	delete[] rsc_blend;

	for (int i = 0; i < ORG_C_ROWS; i++) {
		delete[] dfc[i];
		delete[] fc_idata[i];
		delete[] fc_pdata[i];
		delete[] dcc[i];
		delete[] cc_idata[i];
		delete[] cc_pdata[i];
		delete[] rs[i];
		
	}
	delete[] dfc;
	delete[] fc_idata;
	delete[] fc_pdata;
	delete[] dcc;
	delete[] cc_idata;
	delete[] cc_pdata;
	delete[] rs;

	return 0;
}

// blend_dir	forward or backward
// snesor_p		sensor data pairs
// pimage		prediction images
// pidx			prediction index 0, 1, ..., n predictions

// the following function for block processing is not tested and used yet.

int PSRFM::psrfm_blending_block(int blend_dir, int pidx, SENSOR_PAIR *sensor_p[], short ***pimage, CONTROL_PARAMETER *ctrl, IMGAG_DATE_INFO *image_date_info, short ***blend_output, short ***blend_uncertainty, string fine_sensor_name, ofstream& log_file) {

//	clock_t begin;
//	clock_t end;
//	double elapsed_secs;

	int k;						// for cluster loop
	int b;						// for band loop

	unsigned char **iClusters;
	unsigned char **saved_iClusters;

	short ***saved_blend_image;
	unsigned short ***saved_blend_image_uncertainty;


	ofstream prd_output;			// prediction output stream
	string prd_output_fname;		// prediction output file name

	ofstream unctn_output;			// uncertainty output stream
	string unctn_output_fname;		// uncertainty output file name

	ofstream cluster_output;		// cluster output stream
	string cluster_output_fname;	// cluster output file name

	int ORG_C_ROWS;				// original number of rows of coarse image 
	int ORG_C_COLS;				// original number of columns of coarse image 

	ORG_C_ROWS = ctrl->NUM_ROWS / ctrl->BLOCK_SIZE;
	ORG_C_COLS = ctrl->NUM_COLS / ctrl->BLOCK_SIZE;

	// for quality assessment
	double CC = 0.0;					// correlation at coarse scale of fine difference and coarse difference
	double CF = 0.0;					// correlation at fine scale of fine difference and coarse difference
	double RR = 0.0;

	double sum_CC = 0.0;				// correlation at coarse scale of fine difference and coarse difference
	double sum_CF = 0.0;				// correlation at fine scale of fine difference and coarse difference
	double sum_RR = 0.0;

	double cluster_optimal_value;
	double cluster_optimal_value_old;


	int max_k = 0;

	float **fc_idata;
	float **fc_pdata;

	int tmp_count = 1;
	double RS2;
	double RS2_ave = 0.0;

	float **dcc;				// different coarse images at coarse resolution
	float **dfc;				// different fine image at coarse resolution
	float **rs;

	float **dff;				// different fine at fine resolution
	float **dcf;				// different coarse at fine resolution

	dff = new float *[ctrl->NUM_ROWS];
	dcf = new float *[ctrl->NUM_ROWS];
	for (int i = 0; i < ctrl->NUM_ROWS; i++) {
		dff[i] = new float[ctrl->NUM_COLS];
		dcf[i] = new float[ctrl->NUM_COLS];
	}

	fc_idata = new float *[ORG_C_ROWS];
	fc_pdata = new float *[ORG_C_ROWS];
	dcc = new float *[ORG_C_ROWS];
	dfc = new float *[ORG_C_ROWS];
	rs = new float *[ORG_C_ROWS];
	for (int i = 0; i < ORG_C_ROWS; i++) {
		fc_idata[i] = new float[ORG_C_COLS];
		fc_pdata[i] = new float[ORG_C_COLS];
		dcc[i] = new float[ORG_C_COLS];
		dfc[i] = new float[ORG_C_COLS];
		rs[i] = new float[ORG_C_COLS];
	}

	saved_blend_image = new short **[ctrl->NUM_BANDS];
	saved_blend_image_uncertainty = new unsigned short **[ctrl->NUM_BANDS];
	for (b = 0; b < ctrl->NUM_BANDS; b++) {
		saved_blend_image[b] = new short *[ctrl->NUM_ROWS];
		saved_blend_image_uncertainty[b] = new unsigned short *[ctrl->NUM_ROWS];
		for (int i = 0; i < ctrl->NUM_ROWS; i++) {
			saved_blend_image[b][i] = new short[ctrl->NUM_COLS];
			saved_blend_image_uncertainty[b][i] = new unsigned short[ctrl->NUM_COLS];
		}
	}

	saved_iClusters = new unsigned char *[ctrl->NUM_ROWS];
	for (int i = 0; i < ctrl->NUM_ROWS; i++) {
		saved_iClusters[i] = new unsigned char[ctrl->NUM_COLS];
	}

	int mk;				// total number of coarse resolution image pixels;

	int BLOCK_ROWS = ctrl->BLOCK_SIZE;
	int BLOCK_COLS = ctrl->BLOCK_SIZE;

	float dt;			// temporal interval of input and prediction images
	int nbThreads;


	mk = ORG_C_ROWS * ORG_C_COLS;


	short **blended_image;
	blended_image = new short *[ctrl->NUM_ROWS];
	for (int i = 0; i < ctrl->NUM_ROWS; i++) {
		blended_image[i] = new short[ctrl->NUM_COLS];
	}


	float **A;				// cluster fraction of a block (a coarse resolution image pixel)
	float *cc_refl_change;


	cc_refl_change = new float[ORG_C_ROWS * ORG_C_COLS];


	float **cc_idata = new float *[ORG_C_ROWS];		// coarse image in coarse resoluton
	float **cc_pdata = new float *[ORG_C_ROWS];		// coarse prediction image in coarse resolution
	for (int i = 0; i < ORG_C_ROWS; i++) {
		cc_idata[i] = new float[ORG_C_COLS];
		cc_pdata[i] = new float[ORG_C_COLS];
	}


	int im;			// im = 0 forward. im = 1 backward


	// for OpenMP
	// get number of cores and use all but minus 1
	int numCores = omp_get_num_procs();
	numCores--;

	Eigen::initParallel();
	nbThreads = Eigen::nbThreads() - 1;
	setNbThreads(nbThreads);


	// output file name is sensor name of fine image + like + prediction date + fwd (bwd) + .dat
	prd_output_fname = ctrl->OUTPUT_DIR;
	prd_output_fname.append("\\");
	prd_output_fname.append(fine_sensor_name);
	prd_output_fname.append("_like_");
	prd_output_fname.append(image_date_info->pdt_date[pidx]);

	unctn_output_fname = prd_output_fname;

	if (blend_dir == FORWARD) {
		im = 0;
		prd_output_fname.append("_FWD.dat");
		unctn_output_fname.append("_FWD_Q.dat");
	}
	else if (blend_dir == BACKWARD) {
		im = 1;
		prd_output_fname.append("_BWD.dat");
		unctn_output_fname.append("_BWD_Q.dat");
	}

	iClusters = new unsigned char *[ctrl->NUM_ROWS];

	for (int i = 0; i < ctrl->NUM_ROWS; i++) {
		iClusters[i] = new unsigned char[ctrl->NUM_COLS];

	}

	prd_output.open(prd_output_fname.c_str(), ios::out | ios::binary);
	if (!prd_output) {
		cout << "Cannot open file for writing: <" << prd_output_fname.c_str() << ">" << endl;
	}

	unctn_output.open(unctn_output_fname.c_str(), ios::out | ios::binary);
	if (!unctn_output) {
		cout << "Cannot open file for writing: <" << unctn_output_fname.c_str() << ">" << endl;
	}

	k = ctrl->CLUSTER_RANGE[0];		// initial number of clusters

	MatrixXf EigenA = MatrixXf::Zero(mk, k);			// A Ax = b
	VectorXf Eigenlv = VectorXf::Zero(mk);				// b
	MatrixXf EigenA2 = MatrixXf::Zero(mk, k);			// weighted A
	VectorXf Eigenlv2 = VectorXf::Zero(mk);				// weighted lv

	VectorXf Eigenfv = VectorXf::Zero(k);				// x
	MatrixXf EigenPv = MatrixXf::Zero(mk, mk);			// weights

														//MatrixXf EigenNv = MatrixXf::Zero(mk, mk);
	MatrixXf EigenNv = MatrixXf::Zero(k, k);		// k cluster number
	MatrixXf EigenQv = MatrixXf::Zero(mk, mk);
	MatrixXf EigenQvv = MatrixXf::Zero(mk, mk);
	MatrixXf EigenRv0 = MatrixXf::Zero(mk, mk);
	MatrixXf EigenRv = MatrixXf::Zero(mk, mk);
	MatrixXf EigenRvDiag = MatrixXf::Zero(mk, mk);

	VectorXf Eigenvv = VectorXf::Zero(mk);

	float sv2;
	float rv;


	// time interval between the input and the prediction date
	dt = float(abs(image_date_info->date_num[im] - image_date_info->pdt_date_num[pidx]));

	cluster_optimal_value_old = 1000000.0;			// cluster optimal value holder

	// for cluster optimization
	for (k = ctrl->CLUSTER_RANGE[0]; k <= ctrl->CLUSTER_RANGE[1]; k++) {

		cout << "Cluster No: " << k << endl;

		log_file << "Cluster No: " << k << endl;

		sum_CC = 0.0;
		sum_CF = 0.0;
		sum_RR = 0.0;

		EigenA.resize(Eigen::NoChange, k);
		EigenA2.resize(Eigen::NoChange, k);
		Eigenfv.resize(k);
		EigenNv.resize(k, k);

		A = new float *[mk];
		for (int ik = 0; ik < mk; ik++) {
			A[ik] = new float[k];
		}

		// clustering - use all the bands for clustering
		KMeans(sensor_p[im]->fimage.data, iClusters, k+1, ctrl->NUM_BANDS, ctrl->NUM_ROWS, ctrl->NUM_COLS, ctrl->ITERATIONS, ctrl->MAX_DIFF_THRESHOLD, ctrl->INVALID);

		// get the cluster percetage within a coarse resolution pixel -- A is the output 
		mk = area_fraction(iClusters, k, ctrl->NUM_ROWS, ctrl->NUM_COLS, ctrl->BLOCK_SIZE, A);

		// for each band
		for (b = 0; b < ctrl->NUM_BANDS; b++) {

			cout << "B = " << b << endl;

			log_file << "B = " << b << " ";

			block_mean(sensor_p[im]->cimage.data[b], cc_idata, pimage[b], cc_pdata, ctrl->NUM_ROWS, ctrl->NUM_COLS, BLOCK_ROWS, BLOCK_COLS, ctrl->INVALID);

			// derive the change rate 
			int ii = 0;
			for (int i = 0; i < ORG_C_ROWS; i++) {
				for (int j = 0; j < ORG_C_COLS; j++) {
					cc_refl_change[ii] = (cc_pdata[i][j] - cc_idata[i][j]) / dt;
					ii++;
				}
			}

			// for weight based on uncertainty input for linear system equations
			double wu = sqrt(2.0)*ctrl->c_error / dt;
			wu = 1.0 / (wu*wu);
			for (int i = 0; i < mk; i++) {
				EigenPv(i, i) = float(wu);
			}

			for (int i = 0; i < mk; i++) {
				Eigenlv(i) = cc_refl_change[i];		// change rate for Eigen operation
				for (int j = 0; j < k; j++) {
					EigenA(i, j) = A[i][j];
				}
			}

			// weighting with uncertainty
			EigenA2.noalias() = EigenPv * EigenA;
			Eigenlv2.noalias() = EigenPv * Eigenlv;
			// EIGEN method to solve linear systems of equations
			Eigenfv = EigenA2.bdcSvd(ComputeThinU | ComputeThinV).solve(Eigenlv2);
			EigenNv.noalias() = EigenA.transpose() * EigenPv * EigenA;
			EigenQv = EigenNv.inverse();

			// the following Rv calculation is replaced by the following few lines due to computation inefficiency
			//EigenRv.noalias() = EigenPv - EigenPv * EigenA*EigenQv*EigenA.transpose()*EigenPv;

			EigenRv0 = EigenA * EigenQv * EigenA.transpose();

			//EigenRv = EigenPv * EigenRv0 * EigenPv;   // take long time...Eigen seems too slow

			// as EigenPv is diagonal matrix, the following is the same as EigenRv = EigenPv * EigenRv0 * EigenPv; 
			// therefore to replace it

			for (int i = 0; i < mk; i++) {
#pragma omp parallel for shared(EigenRv, EigenRv0, EigenPv) num_threads(numCores)
				for (int j = 0; j < mk; j++) {
					EigenRv(i, j) = EigenRv0(i, j)*EigenPv(i, i)*EigenPv(j, j);
				}
			}

			EigenRv.noalias() = EigenPv - EigenRv;			// Rv

			Eigenvv = EigenA * Eigenfv - Eigenlv;			// vv
			EigenRvDiag = EigenRv.diagonal();
			rv = EigenRvDiag.sum();							// rv

			sv2 = Eigenvv.transpose()*EigenPv*Eigenvv;
			sv2 = sv2 / rv;									// sv2

			EigenQvv = sv2 * EigenQv;						// Qvv


			// get the predicted image and its uncertainty
			for (int i = 0; i < ctrl->NUM_ROWS; i++) {
				for (int j = 0; j < ctrl->NUM_COLS; j++) {
					blend_output[b][i][j] = short(sensor_p[im]->fimage.data[b][i][j] + dt * Eigenfv(iClusters[i][j]));
					blend_uncertainty[b][i][j] = short(ctrl->f_error*ctrl->f_error + dt * dt * EigenQvv(iClusters[i][j], iClusters[i][j]));
				}
			}

			// residual & correlation analysis

			block_mean(sensor_p[im]->fimage.data[b], fc_idata, blend_output[b], fc_pdata, ctrl->NUM_ROWS, ctrl->NUM_COLS, BLOCK_ROWS, BLOCK_COLS, ctrl->INVALID);

			tmp_count = 1;
			RS2_ave = 0.0;
			for (int ii = 0; ii < ORG_C_ROWS; ii++) {
				for (int jj = 0; jj < ORG_C_COLS; jj++) {
					dfc[ii][jj] = fc_idata[ii][jj] - fc_pdata[ii][jj];		// fine image difference (input and the blended) in coarse resolution
					dcc[ii][jj] = cc_idata[ii][jj] - cc_pdata[ii][jj];		// coarse image difference (input and the input for prediction) in coarse resolution
					rs[ii][jj] = dcc[ii][jj] - dfc[ii][jj];					// fine and coarse difference in coarse resolution
					RS2 = rs[ii][jj] * rs[ii][jj];
					RS2_ave += (RS2 - RS2_ave) / tmp_count;
					tmp_count++;
				}
			}

			for (int i = 0; i < ctrl->NUM_ROWS; i++) {
				for (int j = 0; j < ctrl->NUM_COLS; j++) {
					dff[i][j] = float(blend_output[b][i][j] - sensor_p[im]->fimage.data[b][i][j]);
					dcf[i][j] = float(pimage[b][i][j] - sensor_p[im]->cimage.data[b][i][j]);
				}
			}

			// correlation of fine image difference (input and the blended - from predicted fine image) (from predicted fine image) and MODIS image for prediction at coarse resolution
			CC = correlation(dfc, dcc, ORG_C_ROWS, ORG_C_COLS, ctrl->INVALID);
			sum_CC += CC;

			// correlation of fine image and MODIS image for prediction
			CF = correlation(dff, dcf, ctrl->NUM_ROWS, ctrl->NUM_COLS, ctrl->INVALID);
			sum_CF += CF;

			// residues
			RR = sqrt(RS2_ave);
			sum_RR += RR;

			cout << "CC = " << CC << "; CF = " << CF << ", RR = " << RR << endl;

			log_file << "CC = " << CC << "; CF = " << CF << ", RR = " << RR << endl;

		}  // end of band loop

		if (!ctrl->CLUSTER_OPTIMAL.compare("CC")) {
			cluster_optimal_value = 1.0 - sum_CC / ctrl->NUM_BANDS;
		}
		else if (!ctrl->CLUSTER_OPTIMAL.compare("CF")) {
			cluster_optimal_value = 1.0 - sum_CF / ctrl->NUM_BANDS;
		}
		else if (!ctrl->CLUSTER_OPTIMAL.compare("RR")) {
			cluster_optimal_value = sum_RR;
		}

		if (cluster_optimal_value < cluster_optimal_value_old) {
			cluster_optimal_value_old = cluster_optimal_value;

			max_k = k;				// max value with k cluster

			for (int b = 0; b < ctrl->NUM_BANDS; b++) {
				for (int i = 0; i < ctrl->NUM_ROWS; i++) {
					for (int j = 0; j < ctrl->NUM_COLS; j++) {
						saved_blend_image[b][i][j] = blend_output[b][i][j];
						saved_blend_image_uncertainty[b][i][j] = (unsigned short)sqrt(blend_uncertainty[b][i][j]);
					}
				}
			}

			for (int i = 0; i < ctrl->NUM_ROWS; i++) {
				for (int j = 0; j < ctrl->NUM_COLS; j++) {
					saved_iClusters[i][j] = iClusters[i][j];
				}
			}
		}

		for (int i = 0; i < k; i++) {
			delete[] A[i];
		}
		delete[] A;

	} // end of k (clustering) loop

	if (!ctrl->CLUSTER_OPTIMAL.compare("CC")) {
		log_file << "Cluster number: " << max_k << " has the maximun average value: CC = " << sum_CC / ctrl->NUM_BANDS << endl;
	}
	else if (!ctrl->CLUSTER_OPTIMAL.compare("CF")) {
		log_file << "Cluster number: " << max_k << " has the maximun average value: CF = " << sum_CF / ctrl->NUM_BANDS << endl;
	}
	else if (!ctrl->CLUSTER_OPTIMAL.compare("RR")) {
		log_file << "Cluster number: " << max_k << " has the smallest average value: RR = " << sum_RR / ctrl->NUM_BANDS << endl;
	}



	// output the optimized blended image
	for (int b = 0; b < ctrl->NUM_BANDS; b++) {
		for (int i = 0; i < ctrl->NUM_ROWS; i++) {
			prd_output.write(reinterpret_cast<char *>(saved_blend_image[b][i]), sizeof(short)*ctrl->NUM_COLS);
			unctn_output.write(reinterpret_cast<char *>(saved_blend_image_uncertainty[b][i]), sizeof(unsigned short)*ctrl->NUM_COLS);
		}
	}
	prd_output.close();
	unctn_output.close();

	cluster_output_fname = ctrl->OUTPUT_DIR; //    "D:\\cluster_";	// cluster output file name
	cluster_output_fname.append("\\");
	cluster_output_fname.append(fine_sensor_name);
	cluster_output_fname.append("_cluster_");
	//cluster_output_fname.append(image_date_info->pdt_date[pidx]);
	if (blend_dir == FORWARD) {
		cluster_output_fname.append(image_date_info->f_date[0]);
		cluster_output_fname.append("_FWD_");
	}
	else {
		cluster_output_fname.append(image_date_info->f_date[1]);
		cluster_output_fname.append("_BWD_");
	}

	cluster_output_fname.append(to_string(max_k));
	cluster_output_fname.append(".dat");

	cluster_output.open(cluster_output_fname.c_str(), ios::out | ios::binary);

	// output cluster of the best image
	for (int i = 0; i < ctrl->NUM_ROWS; i++) {
		cluster_output.write(reinterpret_cast<char *>(saved_iClusters[i]), sizeof(unsigned char)*ctrl->NUM_COLS);
	}
	cluster_output.close();

	// return the optimized blended images
	for (int b = 0; b < ctrl->NUM_BANDS; b++) {
		for (int i = 0; i < ctrl->NUM_ROWS; i++) {
			for (int j = 0; j < ctrl->NUM_COLS; j++) {
				blend_output[b][i][j] = saved_blend_image[b][i][j];
				blend_uncertainty[b][i][j] = saved_blend_image_uncertainty[b][i][j] * saved_blend_image_uncertainty[b][i][j];
			}
		}
	}

	// free memory

	for (int b = 0; b < ctrl->NUM_BANDS; b++) {
		for (int i = 0; i < ctrl->NUM_ROWS; i++) {
			delete[] saved_blend_image[b][i];
			delete[] saved_blend_image_uncertainty[b][i];
		}
		delete[] saved_blend_image[b];
		delete[] saved_blend_image_uncertainty[b];
	}
	delete[] saved_blend_image;
	delete[] saved_blend_image_uncertainty;


	for (int i = 0; i < ctrl->NUM_ROWS; i++) {
		delete[] iClusters[i];
		delete[] dff[i];
		delete[] dcf[i];
	}

	delete[] iClusters;
	delete[] dff;
	delete[] dcf;

	for (int i = 0; i < ORG_C_ROWS; i++) {
		delete[] dfc[i];
		delete[] fc_idata[i];
		delete[] fc_pdata[i];
		delete[] dcc[i];
		delete[] cc_idata[i];
		delete[] cc_pdata[i];
		delete[] rs[i];

	}
	delete[] dfc;
	delete[] fc_idata;
	delete[] fc_pdata;
	delete[] dcc;
	delete[] cc_idata;
	delete[] cc_pdata;
	delete[] rs;

	return 0;
}


// pidx: prediction index 
// fine_sensor_name: for output file name prefix
int PSRFM::psrfm_merge_output(int pidx, string fine_sensor_name, short ***forward_output, short ***forward_uncertainty, short ***backward_output, short ***backward_uncertainty, CONTROL_PARAMETER *ctrl, IMGAG_DATE_INFO *image_date_info) {

	ofstream blend_output_file;			// blended prediction output stream
	string blend_output_fname;
	string blend_envi_hdr_fname;

	ofstream blend_unctn_file;			// blended uncertainty 
	string blend_unctn_fname;

	int i;
	int j;
	int b;

	short **blend_output;
	unsigned short **blend_uncertainty;

	float dt01;
	float dt21;
	float w01;
	float w21;
	float wsum;
	
	blend_output = new short *[ctrl->NUM_ROWS];
	blend_uncertainty = new unsigned short *[ctrl->NUM_ROWS];
	for (i = 0; i < ctrl->NUM_ROWS; i++) {
		blend_output[i] = new short[ctrl->NUM_COLS];
		blend_uncertainty[i] = new unsigned short[ctrl->NUM_COLS];
	}
	
	blend_output_fname = ctrl->OUTPUT_DIR;
	blend_output_fname.append("\\");
	blend_output_fname.append(fine_sensor_name);
	blend_output_fname.append("_like_");
	blend_output_fname.append(image_date_info->pdt_date[pidx]);

	blend_envi_hdr_fname = blend_output_fname;
	blend_unctn_fname = blend_output_fname;
	blend_unctn_fname.append("_Q.dat");

	blend_output_fname.append("_PSRFM.dat");
	blend_envi_hdr_fname.append("_PSRFM.hdr");

	blend_output_file.open(blend_output_fname.c_str(), ios::out | ios::binary);
	blend_unctn_file.open(blend_unctn_fname.c_str(), ios::out | ios::binary);
	
	if (!ctrl->MERGE_METHOD.compare("temporal")) {   // merge based on time interval
		
		dt01 = float(abs(image_date_info->date_num[0] - image_date_info->pdt_date_num[pidx]));
		dt21 = float(abs(image_date_info->date_num[1] - image_date_info->pdt_date_num[pidx]));

		w01 = float(1.0 / dt01);
		w21 = float(1.0 / dt21);

		wsum = w01 + w21;
		w01 = w01 / wsum;
		w21 = w21 / wsum;

		for (b = 0; b < ctrl->NUM_BANDS; b++) {
			for (i = 0; i < ctrl->NUM_ROWS; i++) {
				for (j = 0; j < ctrl->NUM_COLS; j++) {
					if (forward_output[b][i][j] != ctrl->INVALID && backward_output[b][i][j] != ctrl->INVALID) {
						blend_output[i][j] = short(forward_output[b][i][j] * w01 + backward_output[b][i][j] * w21);
						blend_uncertainty[i][j] = (unsigned short)(forward_uncertainty[b][i][j] * w01 * w01 + backward_uncertainty[b][i][j] * w21 * w21);
						blend_uncertainty[i][j] = (unsigned short)sqrt(blend_uncertainty[i][j]);
					}
					else if (forward_output[b][i][j] != ctrl->INVALID && backward_output[b][i][j] == ctrl->INVALID) {
						blend_output[i][j] = short(forward_output[b][i][j]);
						blend_uncertainty[i][j] = (unsigned short)(forward_uncertainty[b][i][j]);
						blend_uncertainty[i][j] = (unsigned short)sqrt(blend_uncertainty[i][j]);
					}
					else if (forward_output[b][i][j] == ctrl->INVALID && backward_output[b][i][j] != ctrl->INVALID) {
						blend_output[i][j] = short(backward_output[b][i][j]);
						blend_uncertainty[i][j] = (unsigned short)(backward_uncertainty[b][i][j]);
						blend_uncertainty[i][j] = (unsigned short)sqrt(blend_uncertainty[i][j]);
					}
					else {
						blend_output[i][j] = short(ctrl->INVALID);
						blend_uncertainty[i][j] = (unsigned short)(ctrl->INVALID);
					}
				}
				blend_output_file.write(reinterpret_cast<char *>(blend_output[i]), sizeof(short)*ctrl->NUM_COLS);
				blend_unctn_file.write(reinterpret_cast<char *>(blend_uncertainty[i]), sizeof(unsigned short)*ctrl->NUM_COLS);
			}
		}
	}
	else if (!ctrl->MERGE_METHOD.compare("uncertainty")) {		// merge based on uncertainty

		for (b = 0; b < ctrl->NUM_BANDS; b++) {
			for (i = 0; i < ctrl->NUM_ROWS; i++) {
				for (j = 0; j < ctrl->NUM_COLS; j++) {

					w01 = float(1.0 / forward_uncertainty[b][i][j]);
					w21 = float(1.0 / backward_uncertainty[b][i][j]);

					wsum = w01 + w21;
					w01 = w01 / wsum;
					w21 = w21 / wsum;

					if (forward_output[b][i][j] != ctrl->INVALID && backward_output[b][i][j] != ctrl->INVALID) {
						blend_output[i][j] = short(forward_output[b][i][j] * w01 + backward_output[b][i][j] * w21);
						blend_uncertainty[i][j] = (unsigned short)(forward_uncertainty[b][i][j] * w01 * w01 + backward_uncertainty[b][i][j] * w21 * w21);
						blend_uncertainty[i][j] = (unsigned short)sqrt(blend_uncertainty[i][j]);
					}
					else if (forward_output[b][i][j] != ctrl->INVALID && backward_output[b][i][j] == ctrl->INVALID) {
						blend_output[i][j] = short(forward_output[b][i][j]);
						blend_uncertainty[i][j] = (unsigned short)(forward_uncertainty[b][i][j]);
						blend_uncertainty[i][j] = (unsigned short)sqrt(blend_uncertainty[i][j]);
					}
					else if (forward_output[b][i][j] == ctrl->INVALID && backward_output[b][i][j] != ctrl->INVALID) {
						blend_output[i][j] = short(backward_output[b][i][j]);
						blend_uncertainty[i][j] = (unsigned short)(backward_uncertainty[b][i][j]);
						blend_uncertainty[i][j] = (unsigned short)sqrt(blend_uncertainty[i][j]);
					}
					else {
						blend_output[i][j] = short(ctrl->INVALID);
						blend_uncertainty[i][j] = (unsigned short)(ctrl->INVALID);
					}
				}
				blend_output_file.write(reinterpret_cast<char *>(blend_output[i]), sizeof(short)*ctrl->NUM_COLS);
				blend_unctn_file.write(reinterpret_cast<char *>(blend_uncertainty[i]), sizeof(unsigned short)*ctrl->NUM_COLS);
			}
		}
	}

	blend_output_file.close();
	blend_unctn_file.close();

	// copy a ENVI header file for geospatial information
	ifstream source(ctrl->ENVI_HDR, ios::binary);
	ofstream dest(blend_envi_hdr_fname, ios::binary);

	dest << source.rdbuf();

	source.close();
	dest.close();

	for (i = 0; i < ctrl->NUM_ROWS; i++) {
		delete[] blend_output[i];
		delete[] blend_uncertainty[i];
	}

	delete[] blend_output;
	delete[] blend_uncertainty;

	return 0;
}

void PSRFM::doOnePixel(int i, int j, double **cluster_mean, short ***idata, unsigned char **idata_cluster, int K, int B, int inval) {

	double sum;
	double sum_sqrt;
	double dist;
	double dist_min;
	int inval_count = 0;
	int this_cluster = 0;

	dist_min = 100000000.0;

	// check if this pixel is invalid
	/*for (int b = 0; b < B; b++) {
		if (idata[b][i][j] == inval) {
			inval_count++;
		}
	}

	if (inval_count > 0) {
		this_cluster = K - 1;
	}
	else {
		for (int k = 0; k < K; k++) {
			sum = 0.0;
			for (int b = 0; b < B; b++) {
				dist = idata[b][i][j] - cluster_mean[k][b];
				sum += dist * dist;
			}

			sum_sqrt = sqrt(sum);
			if (sum_sqrt < dist_min) {
				dist_min = sum_sqrt;
				this_cluster = k;
			}
		}
	}*/

	for (int k = 0; k < K; k++) {
		sum = 0.0;
		for (int b = 0; b < B; b++) {

			dist = idata[b][i][j] - cluster_mean[k][b];
			sum += dist * dist;
		}

		sum_sqrt = sqrt(sum);
		if (sum_sqrt < dist_min) {
			dist_min = sum_sqrt;
			this_cluster = k;
		}
	}

	idata_cluster[i][j] = this_cluster;
}


void PSRFM::KMeans(short ***idata, unsigned char **idata_cluster, int K, int B, int rows, int cols, int iterations, float MAX_DIFF_THRESHOLD, int inval) {

	//ofstream mean_out;
	//mean_out.open("mean_out.txt", std::ofstream::out | std::ofstream::app);

	int c_id = 0;
	int find_cnt;
	long inval_count = 0;

	int K1;
	int cluster_row;
	int cluster_col;

	int *cluster_rows = new int[K];
	int *cluster_cols = new int[K];

	for (int k = 0; k < K; k++) {
		cluster_rows[k] = 0;
		cluster_cols[k] = 0;
	}

	double **cluster_mean;
	double **cluster_mean_old;

	cluster_mean = new double *[K];
	cluster_mean_old = new double *[K];

	for (int k = 0; k < K; k++) {
		cluster_mean[k] = new double[B];
		cluster_mean_old[k] = new double[B];
	}


	// count grid points with at least one invalid pixel in a band
	inval_count = 0;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if (idata[0][i][j] == inval) {
				inval_count++;
			}
		}
	}
	cout << "Invalid pixels in input data: " << inval_count << endl;

	// deal with inalid pixels and set the last cluster for invalid pixels
	if (inval_count > 1) {
		K1 = K - 1;
		for (int b = 0; b < B; b++) {
			cluster_mean[K - 1][b] = inval;
		}
	}
	else {
		K1 = K;
	}

	// initialize the clusters
	for (int k = 0; k < K1; k++) {
		find_cnt = 0;
		cluster_row = rand() % rows;
		cluster_col = rand() % cols;

		for (int n = 0; n < c_id; n++) {
			if (cluster_row == cluster_rows[n] && cluster_col == cluster_cols[n]) {
				find_cnt = 1;
				break;
			}
		}
		if (find_cnt == 0) {
			cluster_rows[c_id] = cluster_row;
			cluster_cols[c_id] = cluster_col;

			for (int b = 0; b < B; b++) {
				cluster_mean[c_id][b] = idata[b][cluster_row][cluster_col];
			}
			c_id++;
		}
	}

	// save to compare later
	for (int k = 0; k < K; k++) {
		for (int b = 0; b < B; b++) {
			cluster_mean_old[k][b] = cluster_mean[k][b];
		}
	}

	long *cluster_count = new long[K];

	double max_diff = 0.0;
	double center_diff;

	// get number of cores and use all but minus 1
	int numCores = omp_get_num_procs();
	numCores--;


	// assign each pixel to one of the clusters
	for (int i = 0; i < rows; i++) {
#pragma omp parallel for shared(cluster_mean, idata, idata_cluster) num_threads(numCores)
		for (int j = 0; j < cols; j++) {
			if (idata[0][i][j] == inval) {
				idata_cluster[i][j] = K-1;				
			}
			else {
				doOnePixel(i, j, cluster_mean, idata, idata_cluster, K, B, inval);
			}
			// doOnePixel(i, j, cluster_mean, idata, idata_cluster, K, B, inval);
		}
	}
	

	// iteration
	for (int iter = 0; iter < iterations; iter++) {

		//cout << "Iteration = " << iter << endl;
		max_diff = 0;

		for (int k = 0; k < K1; k++) {  // K or K-1 ?
			cluster_count[k] = 0;
			for (int b = 0; b < B; b++) {
				cluster_mean[k][b] = 0.0;
			}
		}

		// calculate new cluster mean for each cluster using the previous cluster information of every pixel
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				c_id = idata_cluster[i][j];
				cluster_count[c_id]++;
				for (int b = 0; b < B; b++) {
					cluster_mean[c_id][b] += (idata[b][i][j] - cluster_mean[c_id][b]) / cluster_count[c_id];
				}
			}
		}

		// calculate the differce after the new iteration

		for (int k = 0; k < K; k++) {
			//cout << "Cluster count " << k << ": " << cluster_count[k] << endl;
			for (int b = 0; b < B; b++) {
				center_diff = abs(cluster_mean_old[k][b] - cluster_mean[k][b]);
				if (center_diff > max_diff) {
					max_diff = center_diff;
				}
			}
		}

		if (max_diff < MAX_DIFF_THRESHOLD) {
			break;
		}

		for (int k = 0; k < K1; k++) {	// K or K-1 ?
			for (int b = 0; b < B; b++) {
				cluster_mean_old[k][b] = cluster_mean[k][b];
			}
		}

		// derive the new cluster id for every pixel using the new cluster mean if difference is big

		for (int i = 0; i < rows; i++) {
#pragma omp parallel for shared(cluster_mean, idata, idata_cluster) num_threads(numCores)
			for (int j = 0; j < cols; j++) {
				if (idata[0][i][j] == inval) {
					idata_cluster[i][j] = K - 1;
				}
				else {
					doOnePixel(i, j, cluster_mean, idata, idata_cluster, K, B, inval);
				}
				// doOnePixel(i, j, cluster_mean, idata, idata_cluster, K, B, inval);
			}
		}

	}  // iteration

	// check count of last cluster for invalid pixels
	inval_count = 0;
	if (K1 == K - 1) {
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				if (idata_cluster[i][j] == K - 1) {
					inval_count++;
				}
			}
		}
	}
	cout << "Invlaid pixels in output cluster " << inval_count << endl;

	delete[] cluster_rows;
	delete[] cluster_cols;
	delete[] cluster_count;

	for (int k = 0; k < K; k++) {
		delete[] cluster_mean[k];
		delete[] cluster_mean_old[k];
	}
	delete[] cluster_mean;
	delete[] cluster_mean_old;

}


void PSRFM::CRatio(short ***idata, unsigned char **idata_cluster, int K, int b, int rows, int cols, int inval) {

	double min = 100000.0;
	double max = -100000.0;
	double delta = 0.0;
	int c_id = 0;
	int K1;
	long inval_cout = 0;

	// find min and max values, reset the invalid values to 0
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if (idata[b][i][j] == inval) {
				idata_cluster[i][j] = K-1;
				inval_cout++;
			}
			if (idata[b][i][j] > max) {
				max = idata[b][i][j];
			}
			if (idata[b][i][j] < min) {
				min = idata[b][i][j];
			}
			idata_cluster[i][j] = 0;
		}
	}
	cout << "Invalid pixels in cluster data " << inval_cout << endl;

	if (inval_cout > 0) {
		K1 = K - 1;
		delta = (max - min) / (K - 1);		// note: the cluster K is for invalid pixels
	}
	else {
		K1 = K;
		delta = (max - min) / K;			// note: the cluster K is for invalid pixels
	}


	// assign cluster numbers
	for (int c = 0; c < K1; c++) {
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				if (idata[b][i][j] >= min+c*delta && idata[b][i][j] < min + (c+1) * delta && idata[b][i][j] != inval) {
					idata_cluster[i][j] = c;
				}
			}
		}
	}

}

void PSRFM::get_image_date_info(string fname, int *date_int, string *date_string) {

	int str_len;
	string year_str;
	string month_str;
	string day_str;

	const int substr_len = 11;

	str_len = int(fname.length());
	int start = str_len - 15;

	*date_string = fname.substr(start, substr_len);

	year_str = (*date_string).substr(7, 4);
	month_str = (*date_string).substr(3, 3);
	day_str = (*date_string).substr(0, 2);

	*date_int = DOY_START[month_str] + stoi(day_str);

	if (stoi(year_str) % 4 == 0) {
		*date_int++;
	}

}

int PSRFM::check_input_parameters(CONTROL_PARAMETER *ctrl, IMGAG_DATE_INFO *image_date_info) {

	if (ctrl->NUM_PREDICTIONS == 0) {
		cout << "No Coarse image input for Prediction" << endl;
		return 1;
	}
	else {
		cout << "Number of Predictions = " << ctrl->NUM_PREDICTIONS << endl;
	}

	for (int i = 0; i < 2; i++) {
		if (image_date_info->f_date[i] != image_date_info->c_date[i]) {
			cout << "Input Fine and Coarse image DATE mismatch at <" << i << ">.." << endl;
			return 1;
		}
	}

	if (!ctrl->CO_REGISTER.compare("YES")) {
		if (ctrl->COARSE_ROWS <= ctrl->NUM_ROWS || ctrl->COARSE_COLS <= ctrl->NUM_COLS) {
			cout << "For co-registration, the dimension of coarse image should bigger than that of fine image." << endl;
			return 1;
		}
	}



	return 0;
}

int PSRFM::area_fraction(unsigned char **clusters, int k, int rows, int cols, int block_size, float **A) {
	
	int mk;

	int i, j;
	int m, n;
	int ik;
	
	int block_count = block_size * block_size;

	int *cluster_count;

	cluster_count = new int[k];
	
	mk = 0;

	for (i = 0; i < rows; i = i + block_size) {
		for (j = 0; j < cols; j = j + block_size) {

			for (ik = 0; ik < k; ik++) {
				cluster_count[ik] = 0;
			}

			for (m = i; m < (i + block_size); m++) {
				for (n = j; n < (j + block_size); n++) {
					cluster_count[clusters[m][n]]++;
				}
			}

			for (ik = 0; ik < k; ik++) {
				A[mk][ik] = float(cluster_count[ik]) / block_count;
			}

			mk++;
		}
	}

	return mk;
}


float PSRFM::bilinear(float q11, float q12, float q21, float q22, int x1, int x2, int y1, int y2, int x, int y) {

	int x2x1, y2y1, x2x, y2y, yy1, xx1;
	x2x1 = x2 - x1;
	y2y1 = y2 - y1;
	x2x = x2 - x;
	y2y = y2 - y;
	yy1 = y - y1;
	xx1 = x - x1;
	return float(1.0 / (x2x1 * y2y1) * (
		q11 * x2x * y2y +
		q21 * xx1 * y2y +
		q12 * x2x * yy1 +
		q22 * xx1 * yy1
		));
}

void PSRFM::residual_adjust(float **rs, float **rsc, int rows, int cols, int block_rows, int block_cols, int inval) {

	int i;
	int j;
	int m;
	int n;

	int x1, x2, y1, y2, x, y;
	float q11, q12, q21, q22;

	for (i = 0; i < rows / block_rows - 1; i++) {
		for (j = 0; j < cols / block_cols - 1; j = j++) {
			q11 = rs[i][j];
			q12 = rs[i][j + 1];
			q21 = rs[i + 1][j];
			q22 = rs[i + 1][j + 1];
			x1 = i * block_rows; x2 = (i + 1)*block_rows;
			y1 = j * block_cols; y2 = (j + 1)*block_cols;

			for (m = x1; m < x2; m++) {
				for (n = y1; n < y2; n++) {
					x = m; y = n;
					if (q11 == inval || q12 == inval || q21 == inval || q22 == inval) {
						rsc[m][n] = 0;
					}
					else {
						if (x == x1 && y == y1) {
							rsc[m][n] = q11;
						}
						else if (x == x2 && y == y1) {
							rsc[m][n] = q21;
						}
						else if (x == x1 && y == y2) {
							rsc[m][n] = q12;
						}
						else if (x == x2 && y == y2) {
							rsc[m][n] = q22;
						}
						else {
							rsc[m][n] = bilinear(q11, q12, q21, q22, x1, x2, y1, y2, x, y);
						}
					}
				}
			}
		}
	}
}


void PSRFM::co_register(short ***ifdata, short ***icdata, short ***ocdata, CONTROL_PARAMETER *ctrl, string out_fname) {

	int dr, dc;
	int i, j, b, m, n, ii, jj;
	int imax, jmax, jjmax, icount;
	float ***fdata;
	float ***cdata;
	double cb;
	double rmax = 0;
	double rjmax = 0;
	int jcount = 0;
	int istart = 0;
	int jstart = 0;

	int B, frows, fcols, inval;

	B = ctrl->NUM_BANDS;
	frows = ctrl->NUM_ROWS;
	fcols = ctrl->NUM_COLS;
	inval = ctrl->INVALID;

	dr = ctrl->COARSE_ROWS - ctrl->NUM_ROWS;
	dc = ctrl->COARSE_COLS - ctrl->NUM_COLS;

	if (dr > 100) {
		istart = dr / 4;
	}
	if (dc > 100) {
		jstart = dr / 4;
	}

	//ofstream coregistered_output;			// prediction output stream
	//string coregistered_output_fname;		// prediction output file name
	//coregistered_output_fname = out_fname;


	// allocate memery for fdata and cdata
	fdata = new float **[B];
	cdata = new float **[B];
	for (b = 0; b < B; b++) {
		fdata[b] = new float *[frows];
		cdata[b] = new float *[frows];
		for (i = 0; i < frows; i++) {
			fdata[b][i] = new float[fcols];
			cdata[b][i] = new float[fcols];
		}
	}

	// copy the fine image data
	for (b = 0; b < B; b++) {
		for (i = 0; i < frows; i++) {
			for (j = 0; j < fcols; j++) {
				fdata[b][i][j] = float(ifdata[b][i][j]);
			}
		}
	}

	// get number of cores and use all but minus 1
	int numCores = omp_get_num_procs();
	numCores--;

	// locate the indces of coarse data which gives the maximun correlation
	rmax = 0; imax = 0; jmax = 0;
	for (i = istart; i < dr; i++) {
		rjmax = 0; jjmax = 0;
		for (j = jstart; j < dc; j++) {
			for (b = 0; b < B; b++) {
				// copy the subset data
				m = 0; 
				for (ii = i; ii < i + frows; ii++) {
					n = 0;
					for (jj = j; jj < j + fcols; jj++) {
						cdata[b][m][n] = float(icdata[b][ii][jj]);
						n++;
					}
					m++;
				}
			}
			cb = 0;
#pragma omp parallel for shared(fdata, cdata) num_threads(numCores)
			for (b = 0; b < B; b++) {
				cb += abs(correlation(fdata[b], cdata[b], frows, fcols, inval));
			}
			cb /= B;
			if (cb > rjmax) {
				rjmax = cb;
				jjmax = j;
				jcount = 0;
			}
			else {
				jcount++;
				if (jcount > 5) {
					break;
				}
			}
		}
		if (rjmax >= rmax) {
			rmax = rjmax;
			imax = i;
			jmax = jjmax;
			icount = 0;
		}
		else {
			icount++;
			if (icount > 5) {
				break;
			}
		}
	}

	// retreive the coarse data for co-registration
	if (ctrl->CO_REGISTER_ROW == 0) {
		ctrl->CO_REGISTER_ROW = imax;
		ctrl->CO_REGISTER_COL = jmax;
		cout << "Co-register row index at start = " << imax << endl;
		cout << "Co-register col index at start = " << jmax << endl;
	}
	else {
		ctrl->CO_REGISTER_ROW += imax;
		ctrl->CO_REGISTER_COL += jmax;
		cout << "Co-register row index at end = " << imax << endl;
		cout << "Co-register col index at end = " << jmax << endl;
		ctrl->CO_REGISTER_ROW /= 2;
		ctrl->CO_REGISTER_COL /= 2;
	}
	cout << "ctrl->CO_REGISTER_ROW = " << ctrl->CO_REGISTER_ROW << endl;
	cout << "ctrl->CO_REGISTER_COL = " << ctrl->CO_REGISTER_COL << endl;

	for (int b = 0; b < B; b++) {
		for (int ii = 0; ii < frows; ii++) {
			for (int jj = 0; jj < fcols; jj++) {
				ocdata[b][ii][jj] = icdata[b][ii + ctrl->CO_REGISTER_ROW][jj + ctrl->CO_REGISTER_COL];
			}
		}
	}

	// save to file
	//coregistered_output.open(coregistered_output_fname.c_str(), ios::out | ios::binary);
	//if (!coregistered_output) {
	//	cout << "Cannot open file for writing: <" << coregistered_output_fname.c_str() << ">" << endl;
	//}
	// output the coregistered image
	//for (int b = 0; b < ctrl->NUM_BANDS; b++) {
	//	for (int i = 0; i < ctrl->NUM_ROWS; i++) {
	//		coregistered_output.write(reinterpret_cast<char *>(ocdata[b][i]), sizeof(short)*ctrl->NUM_COLS);
	//	}
	//}
	//coregistered_output.close();

	// free memory
	for (b = 0; b < B; b++) {
		for (i = 0; i < frows; i++) {
			delete[] fdata[b][i];
			delete[] cdata[b][i];
		}
		delete[] fdata[b];
		delete[] cdata[b];
	}
	delete[] fdata;
	delete[] cdata;
}