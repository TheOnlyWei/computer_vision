// by Wei Shi 
// last modified 09/23/16
// Program 3

#include <cstdio>
#include <iostream>
#include <string>
#include "image.h"

using namespace std;
using namespace ComputerVisionProjects;


//  width = 5 * sigma
int main(int argc, char **argv){

  if (argc!=4) {
    printf("Usage: %s <input image> <output Hough image> <output vote array>\n", argv[0]);
    return 0;
  }

  const string input_file(argv[1]);
  const string output_file(argv[2]);
  const string hough_array_file(argv[3]);


  Image image;

  if (!ReadImage(input_file, &image)) {
    cout <<"Can't open file " << input_file << endl;
    return 0;
  }

	cout << "\n////////////////////////////IMPORTANT!//////////////////////////////\n\n";
	cout << "This program takes as input a binary edge image and outputs its Hough\n";
	cout << "space image and a Hough voting array contained in " << hough_array_file << ".\n";
	cout << "The grayscale values in hough image is relative to the highest\n";
	cout << "number of votes in the hough array.\n\n";
	cout << "hough array vote incrementation: 1\n";
	cout << "rho resolution: 700\n";
	cout << "theta resolution: 700\n\n";
	cout << "////////////////////////////////////////////////////////////////////\n\n";

	Image hough_image;

	// divides pi into theta_resolution pixels
	const size_t rho_resolution = 700;
	const size_t theta_resolution = 700;

	getHoughSpace(&image, &hough_image, 
											  theta_resolution, 
											  rho_resolution,
											  hough_array_file);


  if (!WriteImage(output_file, hough_image)){
    cout << "Can't write to file " << output_file << endl;
    return 0;
  }




	return 0; // main is the only function where omitting return is allowed
}
