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

	CONTROL_PARAMETER *ctrl = new CONTROL_PARAMETER;

	IMGAG_DATE_INFO *image_date_info = new IMGAG_DATE_INFO;

	SENSOR *pimage = new SENSOR[MAX_PREDICTIONS];			//coarse images for prediction 
	string pdata_fname[MAX_PREDICTIONS];

	short ***pdata;					// prediction coarse image data
	short ***forward_output;
	short ***forward_uncertainty;
	short ***backward_output;
	short ***backward_uncertainty;

	string fine_sensor_name;

	PSRFM iPSRFM;

	iPSRFM.initialize_ctrl_parameters(ctrl);

	iPSRFM.initialize_sensors(sensor_p, NPAIRS);

	//iPSRFM.get_inputs(pfname, sensor_p, pimage, ctrl, NPAIRS);
	iPSRFM.get_inputs(pfname, sensor_p, pdata_fname, ctrl, image_date_info, NPAIRS);

	// check input parameters
	iPSRFM.check_input_parameters(ctrl, image_date_info);

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

	// log input data file names
	cout << endl << "Fine image on start date: " << sensor_p[0]->fimage.fname << endl;
	log_file << endl << "Fine image on start date: " << sensor_p[0]->fimage.fname << endl;
	cout << "Coarse image on start date: " << sensor_p[0]->cimage.fname << endl;
	log_file << "Coarse image on start date: " << sensor_p[0]->cimage.fname << endl;

	cout << "Fine image on end date: " << sensor_p[1]->fimage.fname << endl;
	log_file << "Fine image on end date: " << sensor_p[1]->fimage.fname << endl;
	cout << "Coarse image on end date: " << sensor_p[1]->cimage.fname << endl;
	log_file << "Coarse image on end date: " << sensor_p[1]->cimage.fname << endl << endl;

	// read image data (input pairs)
	iPSRFM.read_image_data(sensor_p, ctrl, NPAIRS, log_file);

	// prediction
	for (int pi = 0; pi < ctrl->NUM_PREDICTIONS; pi++) {

		iPSRFM.read_pimage_data(pdata, pdata_fname[pi], ctrl);

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

		// log input data file names
		cout << endl << "Coarse image on prediction date: " << sensor_p[1]->fimage.fname << endl;

		cout << "Found total invalid pixels in coarse image on prediction date: " << inval_count << endl;
		log_file << "Found total invalid pixels in coarse image on prediction date: " << inval_count << endl;
		
		//forward prediction
		cout << endl << "Forward blending....." << endl;
		log_file << endl << "Forward blending....." << endl;
		iPSRFM.psrfm_blending(FORWARD, pi, sensor_p, pdata, ctrl, image_date_info, forward_output, forward_uncertainty, fine_sensor_name, log_file);

		//backword prediction
		cout << endl << "Backward blending....." << endl;
		log_file << endl << "Backward blending....." << endl;
		iPSRFM.psrfm_blending(BACKWARD, pi, sensor_p, pdata, ctrl, image_date_info, backward_output, backward_uncertainty, fine_sensor_name, log_file);

		//merge forward and backward predictions
		iPSRFM.psrfm_merge_output(pi, fine_sensor_name, forward_output, forward_uncertainty, backward_output, backward_uncertainty, ctrl, image_date_info);

	}

	// free memory
	iPSRFM.clean_memory(sensor_p, NPAIRS, ctrl, pimage);

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

	end = clock();

	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	cout << endl << "Total time spent = " << elapsed_secs << endl;

	log_file << endl << "Total time spent = " << elapsed_secs << endl << endl;

	log_file.close();

    return 0;
}

