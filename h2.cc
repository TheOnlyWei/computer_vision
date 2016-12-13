// by Wei Shi 
// last modified 09/23/16
// Program 2

#include <cstdio>
#include <iostream>
#include <string>
#include "image.h"

using namespace std;
using namespace ComputerVisionProjects;


//  width = 5 * sigma
int main(int argc, char **argv){

  if (argc!=4) {
    printf("Usage: %s <input image> <gray-level threshold> <output image>\n", argv[0]);
    return 0;
  }

  const string input_file(argv[1]);
  const string threshold_str(argv[2]);
  const string output_file(argv[3]);

	unsigned short threshold = stoi(threshold_str);

  Image image;

  if (!ReadImage(input_file, &image)) {
    cout <<"Can't open file " << input_file << endl;
    return 0;
  }

	cout << "\n////////////////////////////IMPORTANT!//////////////////////////////\n\n";
	cout << "This program takes as input a strength image and turns it into a\n";
	cout << "binary image using threshold ([0, 255]) of " << threshold << ".\n\n";
	cout << "Recommended binary thresholds ([0, 255]):\n";
	cout << "hough_complex_1.pgm: 30.\n";
	cout << "hough_simple_1.pgm: 53.\n";
	cout << "hough_simple_2.pgm: 48.\n\n";
	cout << "////////////////////////////////////////////////////////////////////\n\n";

	toBinary(threshold, &image);

  if (!WriteImage(output_file, image)){
    cout << "Can't write to file " << output_file << endl;
    return 0;
  }

	return 0; // main is the only function where omitting return is allowed
}
