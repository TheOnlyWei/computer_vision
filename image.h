// By Wei Shi (based on code by Ioannis Stamos)
// last modified 11/10/16
// Class for representing a 2D gray-scale image,
// with support for reading/writing pgm images.
// To be used in Computer Vision class.

#ifndef COMPUTER_VISION_IMAGE_H_
#define COMPUTER_VISION_IMAGE_H_

#include <cmath>
#include <cstdlib>
#include <utility>
#include <iomanip>
#include <string>
#include <set>
#include <list>
#include <fstream>
#include "disjoint_sets.h"

namespace ComputerVisionProjects {

// Class for representing a gray-scale image.
// Sample usage:
//   Image one_image;
//   one_image.AllocateSpaceAndSetSize(100, 200);
//   one_image.SetNumberGrayLevels(255);
//   // Creates an image such that each pixel is 150.
//   for (int i = 0; i < 100; ++i)
//     for (int j = 0; j < 200; ++j)
//       one_image.SetPixel(i, j, 150);
//   WriteImage("output_file.pgm", an_image);
//   // See image_demo.cc for read/write image.
class Image {
 public:
  Image(): num_rows_{0}, num_columns_{0},
	   num_gray_levels_{0}, pixels_{nullptr}{ }

	// initializes empty image
  // Image(const size_t& width);
  Image(const Image &an_image);
  Image& operator=(const Image &an_image);

  ~Image();

  // Sets the size of the image to the given
  // height (num_rows) and columns (num_columns).
  void AllocateSpaceAndSetSize(size_t num_rows, size_t num_columns);

  size_t num_rows() const { return num_rows_; }
  size_t num_columns() const { return num_columns_; }
  size_t num_gray_levels() const { return num_gray_levels_; }

	void allocateFlatGrayScale(size_t num_rows, size_t num_columns, unsigned short gray_level);

  void SetNumberGrayLevels(size_t gray_levels) {
    num_gray_levels_ = gray_levels;
  }

  // Sets the pixel in the image at row i and column j
  // to a particular gray_level.
  void SetPixel(size_t i, size_t j, int gray_level) {
    if (i >= num_rows_ || j >= num_columns_) abort();
      pixels_[i][j] = gray_level;
  }

  unsigned short int GetPixel(size_t i, size_t j) const {
    if (i >= num_rows_ || j >= num_columns_) abort();
      return pixels_[i][j];
  }

	//void update();
	void printGrayscales();

 private:
  void DeallocateSpace();

	DisjointSets labels;
  size_t num_rows_ = 0;
  size_t num_columns_ = 0;
  size_t num_gray_levels_ = 0;
  int **pixels_ = nullptr; // 2D array of pixel shade values
	//bool need_update = false;

	//unordered_map<unsigned short, double> centerX;
	//unordered_map<unsigned short, double> centerY;

	// Requires the objects to be separated by shade
	// Returns moment about x-axis or row
	//void setCenterX();

	// Requires the objects to be separated by shade
	// Returns moment about y-axis or row
	//void setCenterY();
};
// Requires the objects to be separated by shade
// Returns moment about x-axis or row
unordered_map<unsigned short, double> getCenterX(Image* image);
// Requires the objects to be separated by shade
// Returns moment about x-axis or row
unordered_map<unsigned short, double> getCenterY(Image* image);

// Reads a pgm image from file input_filename.
// an_image is the resulting image.
// Returns true if  everyhing is OK, false otherwise.
bool ReadImage(const std::string &input_filename, Image *an_image);

// Writes image an_image into the pgm file output_filename.
// Returns true if  everyhing is OK, false otherwise.
bool WriteImage(const std::string &output_filename, const Image &an_image);

// Draws a line of given gray-level color from (x0,y0) to (x1,y1);
// an_image is the output_image.
// IMPORTANT: (x0,y0) and (x1,y1) can lie outside the image
// boundaries, so SetPixel() should check the coordinates passed to it.
void DrawLine(int x0, int y0, int x1, int y1, int color, Image *an_image);

// converts output image to binary based on threshold
// @param threshold if pixel is greater than threshold, pixel is set to 255
// if pixel is less than threshold, pixel is set to 0
void toBinary(Image* output, unsigned short threshold);


// generates a map from label to grayscale from a disjoint set data structure
// of labels created using the label() method.
// @param labels is the disjoint set of labels
// @param base is the initial grayscale value
// @param increment is the increment to the grayscale value for each new object
unordered_map<size_t, unsigned short> generateShadeMap(const DisjointSets& labels,
																													 unsigned short base,
																													 short increment);

// assumes background shade is 0, but it really doesn't matter
// we don't care about the gray scale setting of each object
// we only care that we want to assign each object a unique value
// @param min is the background shade
// @param max is the grayscale that differentiates an object from background of 0
// @param label_size is the size of the label wanted (use large number for noisy images)
// otherwise you get a segmentation fault
DisjointSets label(Image* output, unsigned short min, unsigned short max,
																											size_t label_size);

// Precondition: image must be a binary image
// @param min is the background's grayscale
// @param max are the objects' grayscale
// @param label_size is the amount of labels requires
// improper label_size leads to segmentation faults
void separateObjectsByShade(Image* image, unsigned short base,
																					short increment,
																					unsigned short min = 0,
																					unsigned short max = 255,
																					size_t label_size = 500);

// Requires image objects to be separated by shade to work properly.
// Returns a mapping from object label to area of that object
// Post-condition: returns HOMOGENOUS area
unordered_map<unsigned short, unsigned long> getArea(Image* image);

// Post-condition: returns NONHOMOGENOUS area
unordered_map<unsigned short, unsigned long> getNonhomogenousArea(Image* labeled_image,
																																	Image* nonhomogenous_image);


unordered_map<unsigned short, vector<double>> getSecondMoments(Image* input);

// Post-condition: returns the angle between the orthogonal to the line
// passing through the center of the object and the positive x-axis
unordered_map<unsigned short, double> getTheta(Image* input);
unordered_map<unsigned short, double> getEmin(Image* input);
unordered_map<unsigned short, double> getEmax(Image* input);
unordered_map<unsigned short, double> getRoundness(Image* input);

// Post-condition: returns probability of target given mu and sigma in a
// normal distribution
double gaussian(const double& mu, const double& sigma, const double& target);

// Pre-condition: input image objects needs to be separated by  grayscale values
// @param input is the image file with objects separated by different grayscale values
// @param output is the image file to be written to.
void drawMoments(Image* input, Image* output, unsigned short dot_size);

// x'sin(t) - y'cos(t) = 0, where x' = x - x_bar and y' = y_bar
// x_bar = (1/A) * (sum of i* bij)
// Pre-condition: input image objects needs to be separated by  grayscale values
// @param input is the image file with objects separated by different grayscale values
// @param output is the image file to be written to.
void drawOrientation(Image* input, Image* output);

// Pre-condition: input image objects needs to be separated by  grayscale values
// @param input is the image file with objects separated by different grayscale values
// @param output is the image file to be written to
// @param labels is the label of objects you want to draw orientation line of
void drawSelectLine(Image* input, Image* output, vector<unsigned short>labels);

// Pre-condition: input image objects needs to be separated by  grayscale values
// @param input is the image file with objects separated by different grayscale values
// @param output is the image file to be written to
// @param labels is the label of objects you want to draw orientation line of
void drawSelectCenter(Image* input, Image* output, vector<unsigned short>labels, size_t dot_size);


//////////////////////////////////////////////////////////
//
//										EDGE FINDING
//
//////////////////////////////////////////////////////////

// Post-condition: does not normalize resulting matrix
vector<vector<int>> deriveGaussian(const vector<vector<int>>& input, vector<unsigned int> kernel);

// Pre-condition: mask is an integer mask
// Post-condition: row convolution operation on input using mask
void smoothRow(Image* input, Image* output, vector<unsigned int> mask);

// Pre-condition: mask is an integer mask
// Post-condition: column convolution operation on input using mask
void smoothCol(Image* input, Image* output, vector<unsigned int> mask);

// Pre-condition: output needs to be the same size as input
// Post-Condition: returns 1D gaussian mask and smoothes input
vector<unsigned int> getGaussianMask(Image* input, double sigma);

// @param tl is the lower bound parameter for hysteresis thresholding
// @param th is the upper bound parameter for hysteresis thresholding
// @param range tl and th take values in between [0, range]
// To disable HYSTERESIS THRESHOLDING, set (tl, th, range) to
// (0.0, 0.0, 1.0)
Image findEdge(Image* input, double sigma, double tl, double th, double range);

// Pre-condition: filters or convolves an input using kernel
vector<vector<double>> getFilteredArray(Image* input, vector<vector<int>> kernel,
																							const double normalize_constant = 1);


// Post-condition: returns angle 0, 45, 90, 135 closest to parameter 'angle'
unsigned short getClosestAngle(double angle);

// simple floating point number modulo operator function
// follows mathematic convention that remainder is always
// positive. E.g., for fmod_math(a, b), where a < 0, then
// we multiply b by an integer n until n*b < a. The
// absolute value of the difference between n*b and a is
// the result of fmod_math(a, b).
double fmod_math(double a, double b);

// Post-condition: suppresses low edge values in the neightborhood
// of high edge values in the direction of the high edge value's
// gradient.
void nonmax_suppress(vector<vector<double>>& strength,
										 const vector<vector<double>>& orientation);

// TODO: implement hysteresis thresholding
// @param th is upper threshold. Algorithm checks adjacent edges
// along th (orthogonal to gradient) for edge values greater than
// tl. If they exist, they are added to the list, and the list of
// edges is added to the set of edges.
// @param tl is lower threshold mentioned above.
/*
set<list<pair<size_t, size_t>>> hysteresis_threshold(
													vector<vector<double>>& strength,
													const vector<vector<double>>& orientation,
													double tl, double th);
*/


// @param rho_resolution is divided into two parts, half above
// the theta axis and half below. rho = rho_resolution/2 is treated
// as the theta axis.
// @param theta_resolution divides pi (since theta is in [0, pi])
// into theta_resolution number of pixels.
// @param hough_array_file: the file path to store the hough vote array
// Post-condition: creates a hough image file using passed parameter
// 'hough_image' and writes the hough array to file called
// hough_array_file. The first line of hough_array_file contains the
// fields:
// <theta_axis, theta resolution, rho resolution, theta increment, rho increment>
void getHoughSpace(Image* edge_image, Image* hough_image,
																					 const size_t theta_resolution,
																					 const size_t rho_resolution,
																					 const string hough_array_file);

// This function generates only the Hough image given an array file
Image getHoughImage(const string& hough_array_file);

// TODO: current patch voting doesn't work very well, needs revising
void getHoughPatchSpace(Image* edge_image, Image* hough_image,
																					 const size_t theta_resolution,
																					 const size_t rho_resolution,
																					 const string hough_array_file);


unordered_map<unsigned short, double> getNonhomogenousCenterX(Image* labeled_hough,
																															Image* original_hough);
unordered_map<unsigned short, double> getNonhomogenousCenterY(Image* labeled_hough,
																														 	Image* original_hough);

// @param threshold used in converting hough image to a binary image
// from an input hough array from a .txt hough_array_file. The binary image
// sections the hough space into areas of largest vote.
// Post-condition: returns a list of <theta, rho> pairs of corresponding
// hough lines.
// Uses post-processing filtering of redundant lines.
list<pair<double, double>> getLineParameters(const unsigned short threshold,
																						 const string hough_array_file,
																						 const double theta_error,
																						 const double rho_error,
																						 const unsigned short graylevel = 255,
																						 const unsigned short base = 30,
																						 const unsigned short increment = 25);

// Post-condition: filters out redundant lines, e.g., for (theta, rho) pairs, the
// (theta1, rho1) and (theta2, rho2) line parameters are redundant if
// |theta2 - theta1| < theta_inc_multiple_error * theta_increment
// and |rho2 - rho1| < rho_inc_multiple_error * rho_increment
// @param hough image is the image  version of the hough array.
list<pair<double, double>> filterRedundantLines(Image* hough, list<pair<double, double>> theta_rho_pairs,
																								size_t theta_axis, double theta_increment, double rho_increment,
																								unsigned int theta_inc_multiple_error,
																								unsigned int rho_inc_multiple_error);

// Optional method for viewing the sectioned result of hough image
Image getSectionedHough(const unsigned short threshold,
											  const string hough_array_file,
											  const unsigned short graylevel = 255,
											  const unsigned short base = 30,
											  const unsigned short increment = 25);

// Post-condition: draws onto input the lines defined by the parameters in
// theta_rho_pairs
void drawDetectedLines(Image* input, list<pair<double, double>> theta_rho_pairs);
/*
void drawEdgeLines(Image* input, const Image* binary_image, list<pair<double, double>> theta_rho_pairs);
*/

// @param image is used to get the dimensions of the image to determine whether
// the intersection point is out of bounds.
list<pair<double, double>> getIntersections(Image* image,
																						list<pair<double, double>> theta_rho_pairs);

//////////////////////////////////////////////////////////
//
//								PHOTOMETRIC STEREO
//
//////////////////////////////////////////////////////////

// Precondition: image is of a lambertian sphere with frontal lighting.
// Postcondition: returns a vector {x_center, y_center, radius} of the sphere image
// used for calibrating the lambertian surface
// Make sure to select the proper threshold, otherwise the program will not work as
// expected. Threshold 86 is the default and shouldn't be changed for sphere0.pgm.
vector<double> getSphereCenterAndRadius(Image* image, const unsigned short& threshold,
                                                      const unsigned short& sphere_label);

// Precondition: threshold for highlight is predetermined, meaning the client should
// know the binary image threshold to get the highlight location. The client should also
// pass the gray value label of the sphere.
pair<double, double> getHighlight(const Image * image, const unsigned short& threshold,
                                                 const unsigned short& sphere_label);

// @param center is the center of the sphere
// @param target is the target point to calculate the distance from center
// used mainly as an auxiliary function for getLightSourceMatrix
double getDistanceFromSphereCenter(pair<double, double> center, pair<double, double> target);

vector<double> getUnitVector(vector<double> v);
double getVectorLength(vector<double> v);

// two formulas to find (z-zc) of a sphere:
// 1. (z-zc) = sqrt(r^2 - (x-xc)^2 - (y-yc)^2)
// 2. <x-xc, y-yc, rsin(arccos(D/r))>
// Formula 1 is used.
// Precondition: threshold for highlight is predetermined, meaning the client should
// know the binary image threshold to get the highlight location. The client should also
// pass the gray value label of the sphere.
// @params threshold1, threshold2, threshold3 are used in image1, image2, and image3
// respectively to select the brightest area on the image. The user should test
// image1, image2, and image3 until they find a good threshold that isolates the
// highlight of the calibration images.
// The default thresholds used are 200, 219, and 249 respectively for the images
// sphere1.pgm, sphere2.pgm, and sphere3.pgm.
vector<vector<double>> getLightSourceMatrix(const Image* image1,
                                            const Image* image2,
                                            const Image* image3,
                                            const unsigned short& threshold1,
                                            const unsigned short& threshold2,
                                            const unsigned short& threshold3,
                                            double centerX, double centerY, double radius);

// Postcondition: returns the needle map image using image 1. Needle image is
// the 2D projection onto an image plane of 3D vectors normal to the object's surface
// @param threshold is used for thresholding selection of pixel to calculate
// the needle vector for.
Image getNeedleImage(const Image* image1, const Image* image2, const Image* image3,
                     const unsigned short& threshold, const unsigned int& step,
                     const vector<vector<double>> lightsource_matrix);
// Postcondition: returns the albedo image using image 1.
// @param threshold is used for thresholding selection of pixel to calculate
// the albedo for.
Image getAlbedoImage(const Image* image1, const Image* image2, const Image* image3,
                     const unsigned short& threshold,
                     const vector<vector<double>> lightsource_matrix);


// Postcondition: returns the inverse of 3x3 'matrix', which is expected to be
// non-singular
vector<vector<double>> invert3x3Matrix(const vector<vector<double>>& matrix);

// @param 'square' is a square matrix with columns matching the number of rows
// in column matrix 'column'
vector<double> matrixMultiplySquareWithColumn(vector<vector<double>> square,
                                              vector<unsigned short> column);

// Precondition: used to find unit normal vector representing the direction of
// surface normal
vector<double> calculateNormalVector(vector<vector<double>> inverse_source_direction,
                                         vector<unsigned short> intensity);

}  // namespace ComputerVisionProjects





#endif  // COMPUTER_VISION_IMAGE_H_
