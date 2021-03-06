// PSRFM_Main.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <map>

#include <Eigen/Dense>

#include "PSRFM.h"

using namespace std;
using namespace Eigen;

int main(int argc, char **argv)
{

	// current date/time based on current system
	struct tm newtime;
	time_t now = time(0);
	localtime_s(&newtime, &now);

	clock_t begin;
	clock_t end;
	double elapsed_secs;

	string line;
	string word;
	int count = 0;

	ofstream log_file;
	string log_file_name;


	if (argc != 2) {
		cout << "Usage: PSRFM_Main.exe input_parameter_file";
		return 1;
	}

	string pfname = argv[1];

	begin = clock();

	SENSOR_PAIR *sensor_p[NPAIRS];
	for (int i = 0; i < NPAIRS; i++) {
		sensor_p[i] = new SENSOR_PAIR;
	}

	SENSOR_PAIR *observ_p[MAX_PREDICTIONS];
	for (int i = 0; i < MAX_PREDICTIONS; i++) {
		observ_p[i] = new SENSOR_PAIR;
	}

	CONTROL_PARAMETER *ctrl = new CONTROL_PARAMETER;

	IMGAG_DATE_INFO *image_date_info = new IMGAG_DATE_INFO;

	//SENSOR *pimage = new SENSOR[MAX_PREDICTIONS];			//coarse images for prediction 
	//string pdata_fname[MAX_PREDICTIONS];

	short ***pdata;						// coarse image data on prediction date
	short ***forward_output;
	short ***forward_uncertainty;
	short ***backward_output;
	short ***backward_uncertainty;

	string fine_sensor_name;
	string output_fnames;
	string ofname;
	vector<string> vec_fnames;			// create vector to hold output file names

	PSRFM iPSRFM;

	iPSRFM.initialize_ctrl_parameters(ctrl);

	iPSRFM.initialize_sensors(sensor_p, NPAIRS);

	iPSRFM.initialize_sensors(observ_p, MAX_PREDICTIONS);

	//iPSRFM.get_inputs(pfname, sensor_p, pimage, ctrl, NPAIRS);
	iPSRFM.get_inputs(pfname, sensor_p, observ_p, ctrl, image_date_info, NPAIRS);

	// check input parameters
	iPSRFM.check_input_parameters(ctrl, sensor_p, image_date_info);

	// get sensor name of input fine image as part of the blended file name
	fine_sensor_name = iPSRFM.get_short_fine_sensor_name(sensor_p[0]->fimage.fname);


	pdata = new short **[ctrl->NUM_BANDS];
	forward_output = new short **[ctrl->NUM_BANDS];
	forward_uncertainty = new short **[ctrl->NUM_BANDS];
	backward_output = new short **[ctrl->NUM_BANDS];
	backward_uncertainty = new short **[ctrl->NUM_BANDS];
	
	for (int b = 0; b < ctrl->NUM_BANDS; b++) {
		pdata[b] = new short *[ctrl->NUM_ROWS];
		forward_output[b] = new short *[ctrl->NUM_ROWS];
		forward_uncertainty[b] = new short *[ctrl->NUM_ROWS];
		backward_output[b] = new short *[ctrl->NUM_ROWS];
		backward_uncertainty[b] = new short *[ctrl->NUM_ROWS];
		for (int i = 0; i < ctrl->NUM_ROWS; i++) {
			pdata[b][i] = new short[ctrl->NUM_COLS];
			forward_output[b][i] = new short[ctrl->NUM_COLS];
			forward_uncertainty[b][i] = new short[ctrl->NUM_COLS];
			backward_output[b][i] = new short[ctrl->NUM_COLS];
			backward_uncertainty[b][i] = new short[ctrl->NUM_COLS];
		}
	}


	log_file_name = ctrl->OUTPUT_DIR;
	log_file_name.append("\\blend_log.txt");
	
	log_file.open(log_file_name.c_str(), ios::out | ios::app);

	if (!log_file.is_open()) {
		cout << "cannot open file for output log info" << endl;
	}

	log_file << "Year: " << 1900 + newtime.tm_year << " Month: " << 1 + newtime.tm_mon << " Day: " << newtime.tm_mday << endl;
	log_file << "Time: " << newtime.tm_hour << ":" << newtime.tm_min << ":" << newtime.tm_sec << endl;

	// read image data (input pairs)
	iPSRFM.read_image_data(sensor_p, ctrl, NPAIRS, log_file);

	// prediction
	for (int pi = 0; pi < ctrl->NUM_PREDICTIONS; pi++) {

		// apply Kalman Filter algorithm (KFRFM) if it is selected
		if (!ctrl->PREDICT_MODEL.compare("KFRFM") && pi > 0 ) {
			// reset the coarse image data at start, the fine image data were reset already in the fusion data output
			sensor_p[0]->cimage.fname = observ_p[pi - 1]->cimage.fname;
			for (int b = 0; b < ctrl->NUM_BANDS; b++) {
				for (int i = 0; i < ctrl->NUM_ROWS; i++) {
					for (int j = 0; j < ctrl->NUM_COLS; j++) {
						sensor_p[0]->cimage.data[b][i][j] = pdata[b][i][j];
					}

				}
			}

			// reset the start date  
			image_date_info->f_date[0] = image_date_info->pdt_date[pi-1];
			image_date_info->c_date[0] = image_date_info->pdt_date[pi-1];
			image_date_info->date_num[0] = image_date_info->pdt_date_num[pi-1];		
		}

		// read the coarse image data on prediction date
		iPSRFM.read_pimage_data(pdata, observ_p[pi]->cimage.fname, ctrl);
	
		// check and reset invalid pixel values in MODIS image on prediction date
		long inval_count = 0;
		for (int b = 0; b < ctrl->NUM_BANDS; b++) {
			for (int i = 0; i < ctrl->NUM_ROWS; i++) {
				for (int j = 0; j < ctrl->NUM_COLS; j++) {
					if (pdata[b][i][j] < 0 ||
						pdata[b][i][j] > sensor_p[0]->cimage.scale ||
						pdata[b][i][j] == sensor_p[0]->cimage.fillv) {
						pdata[b][i][j] = ctrl->INVALID;
						inval_count++;
					}
				}
			}
		}

		// log input coarse image file names
		cout << "Number of cores used for processing: " << ctrl->NUM_CORES << endl;
		cout << endl << "Coarse image on prediction date: " << observ_p[pi]->cimage.fname << endl;
		log_file << "Coarse image on prediction date: " << observ_p[pi]->cimage.fname << endl;
		cout << "Invalid pixels in coarse image on prediction date: " << inval_count << endl;
		log_file << "Invalid pixels in coarse image on prediction date: " << inval_count << endl << endl;
		
		//forward prediction
		cout << endl << "Forward blending....." << endl;
		log_file << endl << "Forward blending....." << endl;
		iPSRFM.psrfm_blending(FORWARD, pi, sensor_p, pdata, ctrl, image_date_info, forward_output, forward_uncertainty, fine_sensor_name, log_file);

		//backword prediction
		cout << endl << "Backward blending....." << endl;
		log_file << endl << "Backward blending....." << endl;
		iPSRFM.psrfm_blending(BACKWARD, pi, sensor_p, pdata, ctrl, image_date_info, backward_output, backward_uncertainty, fine_sensor_name, log_file);

		//merge forward and backward predictions
		output_fnames = iPSRFM.psrfm_merge_output(pi, sensor_p, observ_p, fine_sensor_name, forward_output, forward_uncertainty, backward_output, backward_uncertainty, ctrl, image_date_info, log_file);

		cout << endl << "Completed blending. " << endl << endl;
		log_file << "Completed blending. " << endl << endl;
	}

	// free memory
	iPSRFM.clean_memory(sensor_p, NPAIRS, ctrl, observ_p);

	for (int b = 0; b < ctrl->NUM_BANDS; b++) {
		for (int i = 0; i < ctrl->NUM_ROWS; i++) {
			delete[] pdata[b][i];
			delete[] forward_output[b][i];
			delete[] forward_uncertainty[b][i];
			delete[] backward_output[b][i];
			delete[] backward_uncertainty[b][i];
		}
		delete[] pdata[b];
		delete[] forward_output[b];
		delete[] forward_uncertainty[b];
		delete[] backward_output[b];
		delete[] backward_uncertainty[b];
	}
	delete[] pdata;
	delete[] forward_output;
	delete[] forward_uncertainty;
	delete[] backward_output;
	delete[] backward_uncertainty;

	delete ctrl;

	end = clock();

	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	cout << endl << "Total time spent = " << elapsed_secs << endl;

	log_file << endl << "Total time spent = " << elapsed_secs << endl << endl;

	log_file.close();

    return 0;
}

