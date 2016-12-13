// by Wei Shi
// last modified 11/10/16
// s2 1

#include <cstdio>
#include <iostream>
#include <string>
#include <fstream>
#include "image.h"

using namespace std;
using namespace ComputerVisionProjects;

vector<double> readFromFile(const string& file_name);
void writeToFile(const string& file_name,
                 const vector<vector<double>>& direction_data);


int main(int argc, char **argv){

  if (argc!=6) {
    printf("Usage: %s <input parameter file> <image1> <image2> <image3> <output directions file>\n", argv[0]);
    return 0;
  }

  const string input_param_file(argv[1]);
  const string image1_file(argv[2]);
  const string image2_file(argv[3]);
  const string image3_file(argv[4]);
  const string output_param_file(argv[5]);

  Image image1, image2, image3;

  if (!ReadImage(image1_file, &image1)) {
    cout <<"Can't open file " <<image1_file << endl;
    return 0;
  }
  if (!ReadImage(image2_file, &image2)) {
    cout <<"Can't open file " << image2_file << endl;
    return 0;
  }
  if (!ReadImage(image3_file, &image3)) {
    cout <<"Can't open file " << image3_file << endl;
    return 0;
  }

	cout << "\n////////////////////////////IMPORTANT!//////////////////////////////\n\n";
	cout << "This program calculates the light source direction of sphere1.pgm, \n";
	cout << "sphere2.pgm, and sphere3.pgm using the following binary threshold\n";
	cout << "Default binary threshold parameters used to find highlights:\n";
	cout << "sphere1.pgm: 200.\n";
	cout << "sphere2.pgm: 249.\n";
	cout << "sphere3.pgm: 219.\n\n";
	cout << "////////////////////////////////////////////////////////////////////\n\n";

  vector<double> xycenter_radius = readFromFile(input_param_file);
  pair<double, double> center(xycenter_radius[0], xycenter_radius[1]);

  vector<vector<double>> lightsource_matrix = getLightSourceMatrix(
                                              &image1, &image2, &image3,
                                              200, 249, 219,
                                              center.first,
                                              center.second,
                                              xycenter_radius[2]);

  writeToFile(output_param_file, lightsource_matrix);

	return 0;
}
vector<double> readFromFile(const string& file_name) {
  vector<double> x_y_radius;
  fstream in;
  in.open(file_name, fstream::in);

  double input;
  while(in >> input) {
    x_y_radius.push_back(input);
  }

  return x_y_radius;

}

void writeToFile(const string& file_name,
                 const vector<vector<double>>& direction_data) {
  fstream out;
  out.open(file_name, fstream::out);
  out.precision(12);
  for(const auto& i: direction_data) {
    for(const auto& j: i) {
      out << j << ' ';

    }
    out << endl;
  }
  out.close();


}
