#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <limits>				// numeric_limits
#include <float.h>				// DBL_MAX

#include <Eigen/Dense>

#define NPAIRS 2
#define MAX_PREDICTIONS 50		// for the period btw the start and the end dates
#define FORWARD 1  				// forward prediction
#define BACKWARD 2				// backward prediction

using namespace Eigen;
using namespace std;

// template functions used to check invalid values
template<typename T>
bool is_infinite(const T &value)
{
	// Since we're a template, it's wise to use std::numeric_limits<T>
	//
	// Note: std::numeric_limits<T>::min() behaves like DBL_MIN, and is the smallest absolute value possible.
	//

	T max_value = std::numeric_limits<T>::max();
	T min_value = -max_value;

	return !(min_value <= value && value <= max_value);
}

template<typename T>
bool is_nan(const T &value)
{
	// True if NAN
	return value != value;
}

template<typename T>
bool is_valid(const T &value)
{
	return !is_infinite(value) && !is_nan(value);
}

template <typename T>
std::vector<T> GetLinearFit(const std::vector<T>& data, const std::vector<T>& xData)
{
	T xSum = 0, ySum = 0, xxSum = 0, xySum = 0, slope, intercept, residual, sigma = 0;

	for (long i = 0; i < data.size(); i++)
	{
		xSum += xData[i];
		ySum += data[i];
		xxSum += xData[i] * xData[i];
		xySum += xData[i] * data[i];
	}
	slope = (data.size() * xySum - xSum * ySum) / (data.size() * xxSum - xSum * xSum);
	intercept = (xxSum*ySum - xySum * xSum) / (data.size()*xxSum - xSum * xSum);

	for (long i = 0; i < data.size(); i++) {
		residual = (xData[i] * slope + intercept - data[i]);
		sigma += (residual*residual);
	}
	sigma = sqrt(sigma / (data.size() - 2));

	std::vector<T> res;
	res.push_back(slope);
	res.push_back(intercept);
	res.push_back(sigma);
	return res;
}

/* data struct for surface reflectance */
typedef struct {
	string date;				/* observation date */
	int date_num;				/* observation day of year */
	string fname;				/* data file name */
	string vfname;				/* variance file name */
	string cfname;				/* cluster file name */
	string mfname;				/* mask file name */
	int  fillv;					/* fill value used */
	int  range[2];				/* sensor data range */
	float res;					/* spatial resolution */
	float scale;				/* scale factor */
	float uncertainty;			/* data accuracy */
	short ***data;				/* sensor data: i=irow; j=icol */
	unsigned short ***std;		/* standard deviation: i=irow; j=icol */
	unsigned char **cluster;	/* clusters */
	unsigned char **mask;		/* mask data: 0 means VALID, >1 INVALID; i=irow; j=icol */
	short int min;				/* minimum value for ipair */
	short int max;				/* maximum value for ipair */
} SENSOR;

/* struct for input image pair (Landsat and MODIS) */
typedef struct {
	SENSOR fimage;				/* for finer resolution data */
	SENSOR cimage;				/* for coarse resolution data */
} SENSOR_PAIR;

/* struct for control parameters */
typedef struct {
	string OUTPUT_DIR;					// directory for final output
	string TEMP_DIR;					// directory for temperary output data
	string ENVI_HDR;					// ENVI header file for the final output

	int NUM_PAIRS;						/* number of input data pair */
	int NUM_BANDS;						// number of bands
	int NUM_ROWS;						// row numbers of fine image
	int NUM_COLS;						// colum numbers of fine image
	int COARSE_ROWS;					// row numbers of coarse image
	int COARSE_COLS;					// colum numbers of coarse image

	int NUM_PREDICTIONS;				/* number of input prediction dates */
	string PREDICT_MODEL;				// prediction model: (PSRFM, KFRFM)

	string CO_REGISTER;					// option for co-registration (YES, NO)
	int CO_REGISTER_ROW = 0;			// row idex of co-registering coarse image
	int CO_REGISTER_COL = 0;			// col idex of co-registering coarse image

	int RESOLUTION;						// resolution of fine image 
	int INVALID = 32764;				// mark value for invalid value pixels  
	int BLOCK_SIZE;						// resolution ratio of coarse and fine images
	float f_error;						// standard deviation of fine image error
	float c_error;						// standard deviation of coarse image error

	string CLUSTER_METHOD;				// currently only kmeans is implemented
	string CLUSTER_DATA;				// input data used for clustering (fine, fine+coarse, coarse ratio)
	int CLUSTER_RANGE[2];				// cluster number range for optimization
	int CLUSTER_FWD;					// final cluster number used for forward prediction
	int CLUSTER_BWD;					// final cluster number used for backward prediction
	int ITERATIONS = 50;				// for kmean max iterations
	float MAX_DIFF_THRESHOLD = 1.0;		// for kmean iterations
	string CLUSTER_OPTIMAL;				// criteria for cluster o[ptimization (CC, CF...)

	string BIAS_METHOD;					// systematic bias removal method (none, apm, iapm)
	string RC_METHOD;					// residual correction method (none, bliear)

	string MERGE_METHOD;				// fwd and bwd merge weighting method (temporal or uncertainty)
	string MERGE_FINE;					// merge fine image on prediction date if available
	int MERGE_FINE_THRESHOLD = 500;		// difference threadhold (reflectance value) for merging fine image pixels on prediction date
	int MERGE_MASK_EXTEND = 3;			// mask extension (# of pixels) for merging fine image pixels on prediction date

	int NUM_CORES = 1;					// number of cores for parallel processing

} CONTROL_PARAMETER;

typedef struct {
	string f_date[2];					// 0 for start; 1 for end
	string c_date[2];
	int date_num[2];
	string pdt_date[MAX_PREDICTIONS];
	int pdt_date_num[MAX_PREDICTIONS];
} IMGAG_DATE_INFO;



class PSRFM
{
public:
	PSRFM();
	~PSRFM();

	map<string, int> DOY_START{ { "Jan", 0 },{ "Feb", 31 },{ "Mar", 59 },{ "Apr", 90 },{ "May", 120 },{ "Jun", 151 },
	{ "Jul", 181 },{ "Aug", 212 },{ "Sep", 243 },{ "Oct", 273 },{ "Nov", 304 },{ "Dec", 334 } };


	void initialize_ctrl_parameters(CONTROL_PARAMETER *ctrl);

	void initialize_sensors(SENSOR_PAIR *sensor_p[], int PAIRS);

	int read_image_data(SENSOR_PAIR *sensor_p[], CONTROL_PARAMETER *ctrl, int PAIRS, ofstream& log_file);
	int read_pimage_data(short ***pdata, string pdata_fname, CONTROL_PARAMETER *ctrl);
	int read_fimage_data(short ***fdata, string fname, unsigned char **mdata, string mfname, CONTROL_PARAMETER *ctrl);

	//int get_inputs(string ifname, SENSOR_PAIR *sensor_p[], SENSOR *pimage, CONTROL_PARAMETER *ctrl, int pairs);

	int get_inputs(string ifname, SENSOR_PAIR *sensor_p[], SENSOR_PAIR *observ_p[], CONTROL_PARAMETER *ctrl, IMGAG_DATE_INFO *image_date_info, int pairs);

	int check_input_parameters(CONTROL_PARAMETER *ctrl, SENSOR_PAIR *sensor_p[], IMGAG_DATE_INFO *image_date_info);


	int psrfm_blending(int blend_dir, int pidx, SENSOR_PAIR *sensor_p[], short ***pdata, CONTROL_PARAMETER *ctrl, IMGAG_DATE_INFO *image_date_info, short ***forward_output, short ***forward_uncertainty, string fine_sensor_name, ofstream& log_file);
	int psrfm_blending_block(int blend_dir, int pidx, SENSOR_PAIR *sensor_p[], short ***pdata, CONTROL_PARAMETER *ctrl, IMGAG_DATE_INFO *image_date_info, short ***forward_output, short ***forward_uncertainty, string fine_sensor_name, ofstream& log_file);
	
	//int psrfm_merge_output(short ***forward_output, short ***forward_uncertainty, short ***backward_output, short ***backward_uncertainty, CONTROL_PARAMETER *ctrl);
	string psrfm_merge_output(int pidx, SENSOR_PAIR *sensor_p[], SENSOR_PAIR *observ_p[], string fine_sensor_name, short ***forward_output, short ***forward_uncertainty, short ***backward_output, short ***backward_uncertainty, CONTROL_PARAMETER *ctrl, IMGAG_DATE_INFO *image_date_info, ofstream& log_file);
	void clean_memory(SENSOR_PAIR *sensor_p[], int PAIRS, CONTROL_PARAMETER *ctrl, SENSOR_PAIR *observ_p[]);
	
	string get_short_fine_sensor_name(string input_str);

private:

	void KMeans(short ***idata, unsigned char **idata_cluster, CONTROL_PARAMETER *ctrl, int K, int B);

	void CRatio(short ***idata, unsigned char **idata_cluster, int K, int b, int rows, int cols, int inval);

	void doOnePixel(int i, int j, double **cluster_mean, short ***idata, unsigned char **idata_cluster, int K, int B, int inval);
	
	void get_image_date_info(string fname, int *date_int, string *date_string);

	void block_mean(short **idata, float **oidata, short **pdata, float **opdata,  int rows, int cols, int block_rows, int block_cols, int inval);

	int build_equations(short **f_idata, float **cf_idata, short **f_pdata, float **cf_pdata, unsigned char **clusters, float **A, float *l, float dt, int k, CONTROL_PARAMETER *ctrl);

	void block_mean_and_cc_change(float **oidata, short **idata, float **opdata, short **pdata, float *cc_change, int rows, int cols, int block_rows, int block_cols);

	double correlation(float **x1, float **x2, int rows, int cols, int inval);
	
	string string_split(string input_str, char delimiter, int position);

	int area_fraction(unsigned char **clusters, int k, int rows, int cols, int block_size, float **A);

	float bilinear(float q11, float q12, float q21, float q22, int x1, int x2, int y1, int y2, int x, int y);

	void residual_adjust(float **rs, float **rsc, int rows, int cols, int block_rows, int block_cols, int inval);

	void co_register(short ***ifdata, short ***icdata, short ***ocdata, CONTROL_PARAMETER *ctrl, string out_fname);
	
};

