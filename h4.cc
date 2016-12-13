// by Wei Shi 
// last modified 10/09/16
// Program 4


#include <cstdio>
#include <iostream>
#include <string>
#include "image.h"

using namespace std;
using namespace ComputerVisionProjects;

int main(int argc, char **argv){

  if (argc != 5) {
    printf("Usage: %s <input image> <input vote array> <input Hough threshold> <output line image>\n", argv[0]);
    return 0;
  }

  const string input_image(argv[1]);
  const string input_vote_array(argv[2]);
  const string threshold_str(argv[3]);
  const string output_image(argv[4]);

	const unsigned short threshold = stoi(threshold_str);

  Image image;

  if (!ReadImage(input_image, &image)) {
    cout <<"Can't open file " << input_image << endl;
    return 0;
  }

	cout << "\n////////////////////////////IMPORTANT!//////////////////////////////\n\n";
	cout << "This program takes as input a grayscale image, an input vote array,\n";
	cout << "and an input threshold value to create a line image of all lines\n";
	cout << "detected using the algorithm. The threshold value is used to create\n";
	cout << "a binary image from the hough image to create sections of the Hough\n";
	cout << "image that define the parameters (rho, theta) for the line. This\n";
	cout << "program includes line filtering to reduce number of lines per edge.\n\n";
	cout << "Recommended Thresholds ([0, 255]):\n";
	cout << "hough_complex_1.pgm: 114.\n";
	cout << "hough_simple_1.pgm: 120.\n";
	cout << "hough_simple_2.pgm: 78.\n\n";
	cout << "////////////////////////////////////////////////////////////////////\n\n";


	Image labeled_image;
	// getLineParameters includes post-processing filtering of lines
	list<pair<double, double>> theta_rho_pairs = getLineParameters(threshold, input_vote_array, 
																																 8, 4, 255, 1, 1);

	drawDetectedLines(&image, theta_rho_pairs);

	// TODO: Extra credit:
	// find intersection of all lines before drawing the lines
	// then draw the lines based on the points of intersection
	// but how do you know which line to skip?

  if (!WriteImage(output_image, image)){
    cout << "Can't write to file " << output_image << endl;
    return 0;
  }




	return 0; // main is the only function where omitting return is allowed
}
