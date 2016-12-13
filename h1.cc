// by Wei Shi 
// last modified 09/23/16
// Program 1

#include <cstdio>
#include <iostream>
#include <string>
#include "image.h"

using namespace std;
using namespace ComputerVisionProjects;

//  width = 5 * sigma
int main(int argc, char **argv){

  if (argc!=3) {
    printf("Usage: %s <input image> <output image>\n", argv[0]);
    return 0;
  }

  const string input_file(argv[1]);
  const string output_file(argv[2]);

  Image image;

  if (!ReadImage(input_file, &image)) {
    cout <<"Can't open file " << input_file << endl;
    return 0;
  }

	cout << "\n////////////////////////////IMPORTANT!//////////////////////////////\n\n";
	cout << "This program uses CANNY ENHANCER algorithm with NONMAX SUPPRESSION\n";
	cout << "to detect edges and create an output strength edge image.\n";
	cout << "Gaussian smoothing standard deviation used by this program is 1.\n\n";
	cout << "////////////////////////////////////////////////////////////////////\n\n";
  Image image_out;
	Image dx_image, dy_image;

	// The last three parameters are hysteresis threshold values
	// to ignore hysteresis thresholding, enter 0.0, 0.0, 1.0 for
	// the last three parameters. For further information, see
	// the findEdge function in image.h file
	image_out = findEdge(&image, 1, 0.0, 0.0, 1.0);

  if (!WriteImage(output_file, image_out)){
    cout << "Can't write to file " << output_file << endl;
    return 0;
  }

	return 0; // main is the only function where omitting return is allowed
}
