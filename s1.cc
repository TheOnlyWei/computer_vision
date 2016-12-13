// by Wei Shi
// last modified 11/10/16
// s11

#include <cstdio>
#include <iostream>
#include <string>
#include <fstream>
#include "image.h"

using namespace std;
using namespace ComputerVisionProjects;

void writeToFile(const string& file_name, const vector<double>& data);

int main(int argc, char **argv){

  if (argc!=4) {
    printf("Usage: %s <input original image> <input threshold value> <output parameters file>\n", argv[0]);
    return 0;
  }

  const string input_file(argv[1]);
  const string threshold_str(argv[2]);
  const string output_file(argv[3]);

	unsigned int threshold = stoi(threshold_str);
  Image image;

  if (!ReadImage(input_file, &image)) {
    cout <<"Can't open file " << input_file << endl;
    return 0;
  }

	cout << "\n////////////////////////////IMPORTANT!//////////////////////////////\n\n";
	cout << "This program calculates the centroid and radius of the sphere\n";
	cout << "in file \"sphere0.pgm\".\n\n";
	cout << "Make sure to use a good threshold value, otherwise the program won't\n";
	cout << "produce accurate results. Binary threshold value 86 is suggested.\n\n";
	cout << "////////////////////////////////////////////////////////////////////\n\n";


	vector<double> data = getSphereCenterAndRadius(&image, threshold, 50);
  writeToFile(output_file, data);

	return 0; // main is the only function where omitting return is allowed
}

void writeToFile(const string& file_name, const vector<double>& data) {
  fstream out;
  out.open(file_name, fstream::out);
  out.precision(12);
  for(const auto& i: data) {
    out << i << ' ';
  }
  out.close();


}
